/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "aom/aom_encoder.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/binary_codes_writer.h"
#include "aom_dsp/bitwriter_buffer.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/bitops.h"
#include "aom_ports/mem_ops.h"
#include "aom_ports/system_state.h"
#if CONFIG_BITSTREAM_DEBUG
#include "aom_util/debug_util.h"
#endif  // CONFIG_BITSTREAM_DEBUG

#include "av1/common/cdef.h"
#include "av1/common/cfl.h"
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/entropymv.h"
#include "av1/common/mvref_common.h"
#include "av1/common/pred_common.h"
#include "av1/common/reconinter.h"
#include "av1/common/reconintra.h"
#include "av1/common/seg_common.h"
#include "av1/common/tile_common.h"

#include "av1/encoder/bitstream.h"
#include "av1/encoder/cost.h"
#include "av1/encoder/encodemv.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/mcomp.h"
#include "av1/encoder/palette.h"
#include "av1/encoder/segmentation.h"
#include "av1/encoder/tokenize.h"

#define ENC_MISMATCH_DEBUG 0

int total = 0;  // MC 2019
float bit_length0 = 0.f, bit_length1 = 0.f;

static float av1_intra_mode_layer0_weights[54 * 16] = {
  -1.34472859f, -2.52086377f, -1.21907496f, -2.62442231f, -1.36456108f,
  -1.85076857f, -1.75697982f, -2.70074654f, -2.01186728f, -1.47276330f,
  -1.49210894f, -1.52696955f, -0.47444168f, 0.08867505f,  0.09964295f,
  1.26356518f,  0.01439609f,  -0.81688899f, -0.50793839f, -1.23898613f,
  -0.06767552f, -1.93357253f, -2.16731191f, -1.21160233f, -1.63814270f,
  -2.20185733f, -1.97344089f, -1.98948014f, -0.85417777f, -1.77182496f,
  1.12845099f,  -0.37248057f, 1.56314051f,  -0.99744338f, -0.42270893f,
  -1.00813162f, -1.27447009f, -1.39900613f, -1.27575421f, -1.28639388f,
  -0.47771707f, -0.97837448f, -1.43649125f, -1.93330860f, -1.12258101f,
  -1.65930676f, -1.65564024f, -1.69175899f, -1.06844699f, 0.11572217f,
  -1.39114916f, -1.32224143f, -1.26412702f, -1.08656287f, -0.73283285f,
  -0.51661754f, -3.07390547f, -1.81956255f, -2.49317861f, -2.34918404f,
  -3.11862493f, -3.84495616f, -1.97497737f, -1.74563420f, -1.39492238f,
  -1.82890868f, 1.05744183f,  0.27218592f,  -2.87983131f, 0.85960984f,
  -1.56861508f, -0.31561625f, -0.17744578f, -0.26456827f, -1.68759775f,
  -1.84392929f, -1.85461676f, -1.67833781f, -1.20882034f, -0.23726073f,
  -1.26488149f, -0.33092192f, -0.50538093f, -0.93326998f, 0.44742784f,
  0.13446647f,  -0.85056525f, -0.82760453f, -1.10008419f, -0.12455311f,
  -1.32504094f, -2.01124835f, -1.91583014f, -1.65612626f, -2.60822034f,
  -2.33542585f, -1.11933470f, -1.85556376f, -2.24145532f, -1.73677206f,
  -1.99684584f, -2.08142853f, -1.42412758f, 0.33718309f,  -0.78021812f,
  -0.15221675f, -1.59199429f, -0.86873209f, 0.65333110f,  0.01478232f,
  0.58083361f,  -1.41929781f, 0.07092807f,  -0.65887439f, -0.04127261f,
  -0.14167656f, -0.75929362f, 0.09406663f,  -0.10605180f, -0.63521659f,
  -1.25242639f, 0.07191706f,  -0.48803216f, -0.25908542f, -0.92937392f,
  -0.76943487f, 0.03202001f,  -0.68810582f, -0.04246534f, -0.55197316f,
  -0.70237046f, 0.21076342f,  -0.11062458f, -1.66637504f, -1.00620592f,
  0.12303604f,  -2.00321126f, -0.05222737f, 0.04923328f,  -0.01844162f,
  0.28145379f,  0.41233397f,  -0.79741007f, -0.29028648f, -0.33426377f,
  -0.06650374f, -0.29935929f, -0.34789035f, -1.42418647f, -0.43265161f,
  0.37286478f,  -0.99247521f, -1.03913200f, -0.12697119f, -0.10455442f,
  -1.46723247f, -0.43378970f, 0.38922939f,  -0.86732292f, -0.60063630f,
  -1.61712241f, -1.69657886f, -0.50130206f, 0.27482462f,  -1.00946200f,
  -0.57687825f, -1.45710349f, -1.39979291f, -1.35991192f, -1.24734640f,
  0.28397655f,  -0.62428916f, -0.30639529f, -0.35062689f, 0.90501994f,
  0.05572261f,  -1.71796548f, -0.10479453f, -1.35078228f, -0.82040793f,
  -0.07072417f, -0.48470831f, -0.54324561f, -0.43556020f, -0.57757121f,
  -1.40337694f, -1.25641024f, -0.79666072f, 0.11310165f,  -0.41256809f,
  -0.47851104f, -0.61234587f, -0.23643969f, 0.29792047f,  -0.00426866f,
  -0.47860229f, -0.36118966f, -0.38841838f, -0.79471290f, -0.72751516f,
  -0.56346869f, -0.76650256f, -1.30257404f, -0.60020411f, -0.71272135f,
  -0.53190637f, -0.23315226f, -0.41795951f, -0.98146957f, -0.52313983f,
  -0.26767144f, -0.12230621f, -0.94776309f, -0.49388161f, -1.13231790f,
  -0.92857045f, -0.15959336f, 0.09274069f,  -0.15788731f, -0.13210642f,
  0.04486415f,  0.04011234f,  -0.27483502f, 0.06583183f,  -0.29369175f,
  0.20024052f,  0.08689474f,  -0.04064971f, -0.15069926f, -0.13407582f,
  0.04207417f,  -0.03537172f, -0.29742697f, -0.25155050f, -0.10513961f,
  -0.29204109f, -0.32159388f, 0.11779372f,  -0.22457708f, -0.09925192f,
  -0.26645407f, 0.04704744f,  0.01768643f,  -0.01814654f, -0.20310701f,
  0.12903525f,  -0.19134848f, -0.06448050f, -0.18500680f, -0.29734027f,
  -0.08669252f, -0.37062702f, -0.20113024f, 0.12800379f,  0.06264347f,
  0.08655868f,  0.00908094f,  0.25891271f,  0.18672709f,  0.11581429f,
  0.07032672f,  0.02287451f,  0.15500040f,  -0.02873144f, -0.11420770f,
  -0.18702960f, -0.34368554f, -0.39120013f, -0.35553777f, -0.11240676f,
  -0.08334691f, -1.00772655f, 1.30852950f,  -0.50752068f, -0.25562558f,
  -1.16672623f, 0.99889696f,  0.70206785f,  -0.69752890f, -0.16687194f,
  -0.39348686f, 0.34242854f,  -0.49750787f, 0.10074376f,  1.51656866f,
  -0.23953331f, 0.73290938f,  -0.74924314f, 1.23292553f,  -0.68014765f,
  2.89756489f,  0.76971978f,  0.30499285f,  0.36711833f,  1.38558459f,
  2.20140719f,  0.27596843f,  1.19344234f,  0.38183266f,  1.73299563f,
  0.97819823f,  0.11554833f,  1.19975066f,  -0.45347807f, 0.38320065f,
  -1.05009985f, -0.68232358f, -0.80626684f, -0.12446785f, -0.82795584f,
  -1.07732999f, -1.31500149f, -0.39848766f, -0.71871889f, -0.65616035f,
  -0.73908311f, -0.48298442f, -0.85509539f, -0.54139900f, 0.07094471f,
  -1.08199406f, -0.93419671f, -0.39356738f, -0.53288907f, 0.07252379f,
  0.13129555f,  -0.25044954f, -0.02207282f, 0.08911113f,  0.11216488f,
  -0.25219271f, 0.12284350f,  0.08695591f,  -0.15098892f, 0.09886478f,
  0.09947824f,  0.07349555f,  -0.01804564f, -0.29991126f, -0.08365519f,
  -0.40522352f, -0.33316395f, 0.15315408f,  -0.28754997f, 0.12707956f,
  0.04991391f,  -0.00448665f, 0.17044903f,  0.10012599f,  -0.03537137f,
  -0.27252764f, 0.06565844f,  -0.07328983f, -0.23592365f, -0.04750932f,
  -0.01009854f, -0.26819202f, -0.10041718f, -0.24925797f, -0.01470207f,
  0.02828817f,  -0.07055964f, 0.09243051f,  -0.21502468f, -0.05895968f,
  -0.15196623f, 0.19910243f,  0.12040682f,  -0.16529885f, 0.19040072f,
  -0.04791486f, 0.06809521f,  0.01571533f,  -0.06237814f, -0.20398186f,
  0.14470118f,  -0.19726026f, -0.28585020f, -0.77040279f, -1.47656918f,
  -0.73407972f, -1.07946002f, -1.52267420f, -1.96672571f, -1.25127006f,
  -1.84090972f, -1.24413979f, -0.77641112f, -0.98540813f, -0.88099372f,
  0.79214984f,  0.02849560f,  -0.29718027f, -1.08180916f, -0.45376605f,
  -1.06111443f, -0.58808219f, -1.38724327f, -1.08868885f, -1.45387828f,
  -1.71027756f, -0.94263488f, -1.43084431f, -1.37832665f, -1.42872727f,
  -0.86514193f, -1.37226248f, -1.57131779f, 0.78047705f,  0.26139262f,
  0.27424359f,  -1.30422628f, -1.13016510f, -1.23262525f, -1.33468628f,
  -1.36236286f, -0.94464678f, -1.53929138f, -1.45604563f, -1.85073888f,
  -1.59765244f, -1.70858347f, -0.67296207f, -1.26271939f, -0.60370177f,
  -1.30262315f, -1.12775874f, -0.02354132f, -1.17035508f, -0.55129266f,
  -1.66859019f, -2.20631456f, -0.50770915f, -0.22569266f, -0.49634519f,
  -0.56306195f, 0.32834321f,  -0.41368416f, -0.57404053f, -0.25020000f,
  -0.93529081f, -0.09142505f, 0.33000740f,  -0.96098888f, 0.73082972f,
  0.71344846f,  -0.76279825f, -0.39585066f, -0.26021990f, -0.70687115f,
  -0.55357850f, -0.58935547f, -0.13484390f, -0.03701411f, -0.52440238f,
  -0.61895049f, -0.12163462f, -0.25045338f, -0.77317226f, -0.46779951f,
  -0.00776350f, -0.11452179f, -0.30457970f, -0.18444744f, -0.38152173f,
  -0.95297587f, -0.44156343f, -0.74075842f, -0.54920059f, -0.41488054f,
  -0.31761304f, -0.32309172f, -0.50588894f, -0.63026828f, -0.58538473f,
  -0.00036020f, -0.03134747f, -0.51907688f, -0.16268441f, -0.10506085f,
  -0.16398014f, -0.00426275f, -0.95439047f, -1.23025608f, -1.11553502f,
  -1.27707827f, 0.20668361f,  2.75095868f,  -2.77786589f, 0.45574117f,
  -0.16950387f, 1.60079658f,  -0.63479453f, -0.70688331f, 1.16952634f,
  -0.18442729f, 0.76257491f,  -0.17277105f, 0.52869684f,  -0.01584324f,
  -0.73074567f, 0.94434673f,  -1.22657275f, 0.81408077f,  0.15346859f,
  1.27181983f,  -1.01982427f, -1.00160575f, 0.35602319f,  0.99563819f,
  -0.05786438f, -0.70866770f, 0.28710264f,  -0.36494386f, 0.48481515f,
  -0.26076955f, -0.32406750f, -0.09276336f, -0.14124414f, 1.21412826f,
  -0.57378030f, 0.95357472f,  -0.21294415f, -0.43012840f, -0.17922872f,
  -0.53241891f, -0.07553056f, 0.21424463f,  -0.64071310f, -0.18031065f,
  -0.32409003f, -0.44842401f, -0.75476760f, -0.05700389f, -0.31735164f,
  -0.07066900f, -0.27913710f, 0.18721832f,  -0.97802764f, -0.02957790f,
  0.19795598f,  -0.57443744f, -0.45028266f, -1.56952322f, -0.09173884f,
  -0.14701849f, 0.19308759f,  -1.14278781f, -1.34431505f, -1.21125805f,
  -0.98475969f, -0.64508533f, -0.83854079f, 0.12046282f,  0.68813384f,
  -0.65353161f, 0.16347751f,  -0.12711960f, -0.68575865f, -0.44496951f,
  -0.57436162f, -1.68765593f, -0.28022033f, -1.69762146f, 0.85090858f,
  -1.38507891f, -0.54518372f, -1.17894363f, -1.04924166f, -0.48039410f,
  -0.81639498f, 0.38145441f,  -0.53304338f, -1.06128860f, -1.31422007f,
  -1.50000072f, -0.45184904f, -1.08361661f, -0.04669328f, -1.25660741f,
  -0.36283621f, -0.54259694f, 0.46375713f,  -0.30610615f, -0.47751313f,
  -0.63315547f, -0.69658530f, -0.54776353f, -0.33925962f, 0.17146656f,
  -0.65777695f, -1.31738305f, -0.67088926f, -1.06627762f, 0.51427329f,
  -1.55594468f, -0.13527888f, -0.42415521f, -1.39122665f, -0.69910347f,
  -1.20098555f, -0.77979815f, -0.81679147f, -0.92310315f, 0.49405172f,
  -0.60847849f, -0.52412665f, -1.19899046f, -1.44216835f, 0.14354254f,
  -1.17344701f, -0.96788162f, -1.02055824f, 0.07554962f,  -0.45610791f,
  -0.30308682f, -0.58650035f, -0.74346405f, 1.06583011f,  -0.11827105f,
  -0.54066199f, -0.87915474f, -0.10246918f, -0.60115743f, -0.54705590f,
  -1.32689762f, -1.13265264f, -1.22430491f, -1.80739963f, -1.42482400f,
  -0.56749761f, 0.60556000f,  0.02873548f,  -1.43660271f, -0.60526472f,
  -0.63982242f, 0.91066074f,  -0.94624561f, -0.62588245f, -0.24445516f,
  -0.04893875f, -0.01936267f, -0.39006403f, 0.12945826f,  -1.42453527f,
  -1.37524915f, -1.32013714f, -1.20278358f, 0.09852993f,  -0.10315577f,
  -0.08065794f, -0.09477025f, -0.00994839f, 0.09520070f,  -0.11385962f,
  -0.08745061f, 0.13500404f,  0.12493370f,  -0.03029579f, 0.17121498f,
  0.03472127f,  -0.07964243f, -0.42111149f, -0.38417643f, -0.20200245f,
  -0.04138506f, -0.09231593f, -0.00623832f, -0.19416910f, 0.12600832f,
  0.00677175f,  0.00019933f,  -0.35941932f, -0.00968444f, -0.16422185f,
  -0.04422733f, 0.16918653f,  -0.19377936f, -0.13898984f, -0.01833220f,
  -0.23003057f, -0.17819738f, -0.23495907f, -0.42173252f, 0.15039498f,
  -0.13094299f, 0.03989442f,  0.03779310f,  0.13587461f,  -0.07663590f,
  -0.22060695f, -0.05867923f, -0.18538265f, -0.13845098f, -0.00619692f,
  0.09337087f,  -0.32394120f, 0.02457848f,  -0.06754370f, -0.02711296f,
  0.04042626f,  -0.43942660f, -0.34323353f, -2.62199759f, -1.83570397f,
  1.32944989f,  -2.31893587f, -2.53645182f, -2.29581738f, -1.89590538f,
  -0.95517159f, -0.88534534f, -1.56841207f, -1.47041833f, -0.62718809f,
  0.38274142f,  0.73662877f,  -2.71332860f, -1.25357664f, -1.87608600f,
  0.66963959f,  -2.05025673f, -0.76931864f, 0.59751183f,  -1.47793198f,
  -1.48561573f, -1.41216469f, -1.05765295f, -1.82081604f, -1.85063028f,
  -1.74727964f, -2.27852201f, -0.73292232f, -0.32425493f, -0.38929012f,
  0.53970754f,  -1.22150409f, -0.69838774f, 0.21952137f,  -0.12319843f,
  -1.53330731f, -1.93789661f, -0.64637381f, -1.61524105f, -1.90739310f,
  -0.91654372f, -1.10830879f, 0.46787357f,  -0.51012194f, -1.22085381f,
  -0.05654236f, 0.21787801f,  -0.55646527f, 1.39780521f,  -0.74258411f,
  -0.50562680f, 0.58642507f,  -0.12347291f, 0.06090624f,  -0.49891523f,
  2.60997629f,  2.17256474f,  1.69830310f,  0.79645658f,  -0.33113375f,
  0.08793503f,  0.22360852f,  0.41350910f,  -0.40948215f, 0.11727334f,
  0.46790299f,  -0.10259704f, 0.02095604f,  -0.42927021f, 1.10943198f,
  0.51955670f,  0.72278011f,  0.34663031f,  2.65269542f,  2.36444712f,
  2.66535425f,  1.25593972f,  -0.06249829f, 1.11259437f,  1.23355770f,
  0.95411611f,  0.14775027f,  0.30185759f,  0.21184132f,  0.49548969f,
  -0.38551453f, -0.51130390f, 0.72222441f,  0.31987080f,  0.46624702f,
  0.47101074f,  1.46061051f,  1.09653759f,  1.08059680f,  0.57475758f,
  0.80132717f,  0.43392190f,  0.76438504f,  0.73639178f,  0.28356266f,
  0.01133403f,  -0.14571588f, -0.07066654f, -0.82332945f, -0.57959372f,
  1.16494417f,  1.00693476f,  0.96404880f,  3.13813162f,  0.64713180f,
  1.40678394f,  1.82910097f,  2.40371013f,  3.35057473f,  1.60572636f,
  1.56177974f,  1.42432034f,  0.12374314f,  -0.30043319f, 0.33148962f,
  -0.60546994f, -1.40367031f, -0.71812409f, -0.18762909f, 0.12894551f,
  -0.57277095f, 1.71871758f,  -0.61199999f, -0.09970672f, -0.48187163f,
  1.29971576f,  1.43933082f,  0.13949873f,  0.09609528f,  -0.11614522f,
  -0.41707113f, -0.20995651f, 0.38771021f,  -0.06047611f, -0.29909387f,
  -0.34229413f, -0.44959259f, -0.31996962f, -0.53914273f, -0.21611452f,
  -0.72331005f, -0.34350699f, 0.12088889f,  -0.08811252f, -0.40167418f,
  -0.37688088f, -0.54099822f, -0.54206634f, -0.54484433f, -0.14436390f,
  0.17844799f,  0.64716429f,  0.81040901f,  0.52024901f,
};

static float av1_intra_mode_layer0_bias[16] = {
  -1.27912414f, -1.19608593f, -0.19420378f, -0.5512073f,
  -0.24851941f, 0.10224217f,  -0.12179026f, -0.61660045f,
  -0.59693199f, 0.1818967f,   -0.93660569f, -1.32011855f,
  -0.17164111f, -1.5121603f,  -0.21590997f, 0.49126673f,
};

static float av1_intra_mode_layer1_weights[16 * 13] = {
  -0.93713713f, -0.72247851f, -0.57512850f, 0.25330722f,  -0.16792443f,
  0.00914478f,  -0.19614846f, -0.17332828f, 0.37497330f,  -0.04883192f,
  0.20458439f,  -0.32846570f, -0.07908146f, -0.48761493f, 0.12397400f,
  0.30885664f,  -1.03973031f, 0.18358073f,  0.51951230f,  -1.01940572f,
  -0.28107375f, -0.45530263f, -0.06533878f, -0.87249565f, -0.73729581f,
  0.89736307f,  0.32965636f,  -1.48205268f, -0.21134375f, -0.66425467f,
  -0.18030094f, 0.06333504f,  0.36394933f,  -2.01219559f, -0.22701602f,
  -0.92664486f, -0.15063332f, 0.80470657f,  -0.37392810f, 0.52831149f,
  -0.68099320f, -0.31497231f, -0.10493826f, -1.36982584f, 0.14514697f,
  -0.53787994f, -0.30737135f, -0.08023137f, -1.69344330f, 0.24123858f,
  -0.09696062f, 0.06578072f,  0.18124302f,  -0.29671329f, 0.12546958f,
  0.25548944f,  -0.30488476f, -0.09978185f, -0.62171650f, -2.51156449f,
  0.18465972f,  0.07399452f,  -0.01837502f, 0.83986384f,  0.14584824f,
  -1.20062542f, -0.17757705f, -0.41853166f, 0.10438190f,  -0.47176892f,
  -0.19204073f, 0.09671723f,  -1.48929286f, -0.09693855f, -0.40407184f,
  -0.79407156f, -0.22919254f, -0.17656256f, 0.98915893f,  0.25288600f,
  -0.61092955f, -1.55243027f, 0.14574085f,  -0.53518820f, 0.00405850f,
  -0.47968265f, -0.48774663f, -0.16461237f, -0.39927962f, 0.59019256f,
  -0.48605397f, -1.84833467f, -0.26481029f, 0.10532492f,  0.73564571f,
  0.25189033f,  -0.23137036f, -1.13673031f, -0.52030319f, 0.08406730f,
  -0.41513166f, 0.39864582f,  -0.16957380f, -0.41774458f, -0.56979561f,
  -0.28237373f, 0.25127643f,  -1.34483385f, 0.04228702f,  -0.67438042f,
  0.62757015f,  0.22157638f,  -0.02209013f, -1.35486436f, -0.13648897f,
  -0.46782866f, 0.07632390f,  0.35197541f,  0.32072407f,  -0.72559828f,
  -0.34066424f, -0.26818177f, -0.21135806f, -1.96569884f, -0.26668221f,
  -0.22678584f, 0.30192733f,  0.65870321f,  -1.57528245f, -1.25249124f,
  -0.01317876f, -0.39135754f, 0.05933694f,  -0.55542612f, 0.07657004f,
  -0.55040556f, -0.59258777f, 0.35454935f,  0.24493609f,  -2.37706327f,
  0.01247345f,  -0.61311400f, 0.13520665f,  0.84179956f,  0.53101110f,
  -1.58692777f, 0.31170231f,  0.37289950f,  -0.15566790f, -0.04442883f,
  0.01785802f,  -0.19399203f, -0.16809368f, 0.02480495f,  -0.75356859f,
  -1.97311914f, -0.09059076f, -0.07336730f, 0.00837672f,  0.32422444f,
  -1.21216309f, -1.25572038f, 0.21622021f,  -1.34141326f, 0.33663988f,
  -0.13823576f, 0.08230273f,  -0.12521884f, -0.66249084f, 0.39310926f,
  0.02396839f,  -0.65418690f, -0.18744481f, -0.00578189f, -0.09468156f,
  0.34678343f,  -0.34314498f, 0.21342035f,  0.17652756f,  0.41512224f,
  0.02563119f,  0.20970821f,  0.21324217f,  -0.97975236f, -0.71767676f,
  -0.09004290f, 0.48560980f,  -2.96616197f, -0.20263413f, -0.53824097f,
  0.07739066f,  0.31717095f,  1.64678121f,  -0.04863302f, -0.75062311f,
  0.37260833f,  -0.41664913f, -0.15529498f, 0.33293751f,  0.14564076f,
  -1.67840254f, 0.21253279f,  0.58225554f,  -1.98118329f, -0.12075335f,
  0.13517460f,  -0.86058992f, -1.03245497f,
};

static float av1_intra_mode_layer1_bias[13] = {
  1.45945013f,  -0.03856213f, 0.03966871f, -1.62548113f, -2.1468184f,
  -2.00411797f, -1.91693163f, -1.4718616f, -1.5770936f,  0.11027168f,
  -1.05009961f, -1.1120646f,  0.41852057f,
};

static float av1_intra_mode_in[54] = { 0 };
static float av1_intra_mode_layer0_out[16] = { 0 };
static float av1_intra_mode_layer0_dout[16] = { 0 };
static float av1_intra_mode_layer0_dW[54 * 16] = { 0 };
static float av1_intra_mode_layer0_db[16] = { 0 };
static float av1_intra_mode_layer1_out[13] = { 0 };
static float av1_intra_mode_layer1_dout[13] = { 0 };
static float av1_intra_mode_layer1_dW[16 * 13] = { 0 };
static float av1_intra_mode_layer1_db[13] = { 0 };

static NN_CONFIG_V2 av1_intra_mode = {
  0,                  // counter (!!Initialize to 0)
  1,                  // num_hidden_layers
  av1_intra_mode_in,  // feature
  {
      // fc layer setting
      {
          // layer 0
          54,                             // num_inputs
          16,                             // num_outputs
          av1_intra_mode_layer0_weights,  // weights
          av1_intra_mode_layer0_bias,     // bias
          RELU,                           // activation
          av1_intra_mode_layer0_out,      // output
          av1_intra_mode_layer0_dout,     // d_out
          av1_intra_mode_layer0_dW,       // dW
          av1_intra_mode_layer0_db,       // db
      },
      {
          16,  // num_inputs (!!same as num_outputs of last layer)
          13,
          av1_intra_mode_layer1_weights,
          av1_intra_mode_layer1_bias,
          NONE,
          av1_intra_mode_layer1_out,
          av1_intra_mode_layer1_dout,
          av1_intra_mode_layer1_dW,
          av1_intra_mode_layer1_db,
      },
  },
  13,                         // num_outputs
  av1_intra_mode_layer1_out,  // logits (!!same as last layer output)
  SOFTMAX_CROSS_ENTROPY,
};

static INLINE void write_uniform(aom_writer *w, int n, int v) {
  const int l = get_unsigned_bits(n);
  const int m = (1 << l) - n;
  if (l == 0) return;
  if (v < m) {
    aom_write_literal(w, v, l - 1);
  } else {
    aom_write_literal(w, m + ((v - m) >> 1), l - 1);
    aom_write_literal(w, (v - m) & 1, 1);
  }
}

static void loop_restoration_write_sb_coeffs(const AV1_COMMON *const cm,
                                             MACROBLOCKD *xd,
                                             const RestorationUnitInfo *rui,
                                             aom_writer *const w, int plane,
                                             FRAME_COUNTS *counts);

static void write_intra_y_mode_kf(FRAME_CONTEXT *frame_ctx,
                                  const MB_MODE_INFO *mi,
                                  const MB_MODE_INFO *above_mi,
                                  const MB_MODE_INFO *left_mi,
                                  const MB_MODE_INFO *aboveleft_mi,
                                  PREDICTION_MODE mode, aom_writer *w) {
  assert(!is_intrabc_block(mi));
  (void)mi;

  // MC 2019
  (void)aboveleft_mi;
  aom_cdf_prob *icdf = get_y_mode_cdf(frame_ctx, above_mi, left_mi);
  float pdf[INTRA_MODES];
  aom_cdf_prob pre = 1 << 15;

  for (int i = 0; i < INTRA_MODES; ++i) {
    pdf[i] = (float)(pre - icdf[i]) / (1 << 15);
    pre = icdf[i];
  }

  bit_length0 -= logf(pdf[mode]);
  ++total;

  PREDICTION_MODE above, left, aboveleft;
  int8_t above_angle, left_angle, al_angle;
  int above_q, left_q, al_q;
  BLOCK_SIZE above_sb, left_sb, al_sb;
  TX_SIZE above_txs, left_txs, al_txs;
  av1_block_mode(above_mi, &above, &above_angle, &above_q, &above_sb,
                 &above_txs);
  av1_block_mode(left_mi, &left, &left_angle, &left_q, &left_sb, &left_txs);
  av1_block_mode(aboveleft_mi, &aboveleft, &al_angle, &al_q, &al_sb, &al_txs);

  float features[54], logits[INTRA_MODES];
  int pt = 0;

  for (int i = 0; i < INTRA_MODES; ++i) {
    features[pt++] = (above == i) ? 1.0f : 0.0f;
  }
  features[pt++] = (float)above_angle / 3.f;
  if (above_sb == 255) {
    features[pt++] = 0.0f;
    features[pt++] = 0.0f;
  } else {
    features[pt++] = mi_size_wide_log2[above_sb] / 5.f;
    features[pt++] = mi_size_high_log2[above_sb] / 5.f;
  }
  if (above_txs == 255) {
    features[pt++] = 0.0f;
    features[pt++] = 0.0f;
  } else {
    features[pt++] = tx_size_wide_log2[above_txs] / 6.f;
    features[pt++] = tx_size_high_log2[above_txs] / 6.f;
  }

  for (int i = 0; i < INTRA_MODES; ++i) {
    features[pt++] = (left == i) ? 1.0f : 0.0f;
  }
  features[pt++] = (float)left_angle / 3.f;
  if (left_sb == 255) {
    features[pt++] = 0.0f;
    features[pt++] = 0.0f;
  } else {
    features[pt++] = mi_size_wide_log2[left_sb] / 5.f;
    features[pt++] = mi_size_high_log2[left_sb] / 5.f;
  }
  if (left_txs == 255) {
    features[pt++] = 0.0f;
    features[pt++] = 0.0f;
  } else {
    features[pt++] = tx_size_wide_log2[left_txs] / 6.f;
    features[pt++] = tx_size_high_log2[left_txs] / 6.f;
  }

  for (int i = 0; i < INTRA_MODES; ++i) {
    features[pt++] = (aboveleft == i) ? 1.0f : 0.0f;
  }
  features[pt++] = (float)al_angle / 3.f;
  if (al_sb == 255) {
    features[pt++] = 0.0f;
    features[pt++] = 0.0f;
  } else {
    features[pt++] = mi_size_wide_log2[al_sb] / 5.f;
    features[pt++] = mi_size_high_log2[al_sb] / 5.f;
  }
  if (al_txs == 255) {
    features[pt++] = 0.0f;
    features[pt++] = 0.0f;
  } else {
    features[pt++] = tx_size_wide_log2[al_txs] / 6.f;
    features[pt++] = tx_size_high_log2[al_txs] / 6.f;
  }

  av1_nn_predict_v2(features, &av1_intra_mode, logits);
  float scores[INTRA_MODES];
  av1_nn_softmax(logits, scores, INTRA_MODES);
  bit_length1 -= logf(scores[mode]);

  aom_write_symbol(w, mode, get_y_mode_cdf(frame_ctx, above_mi, left_mi),
                   INTRA_MODES);

  /*FILE *fp = fopen("intra_mode.txt", "a");
  if (fp) {
    fprintf(fp, "%d | ", mode);
    for (int i = 0; i < INTRA_MODES; ++i) {
      fprintf(fp, "%.5f, ", pdf[i]);
    }
    fprintf(fp, "\n%d | ", mode);
    for (int i = 0; i < INTRA_MODES; ++i) {
      fprintf(fp, "%.5f, ", scores[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
  }*/

  // write data into file (MC 2019)
  /*PREDICTION_MODE above, left, aboveleft;
  int8_t above_angle, left_angle, al_angle;
  int above_q,left_q,al_q;
  BLOCK_SIZE above_sb, left_sb, al_sb;
  TX_SIZE above_txs,left_txs,al_txs;
  av1_block_mode(above_mi, &above, &above_angle,&above_q,&above_sb,&above_txs);
  av1_block_mode(left_mi, &left, &left_angle,&left_q,&left_sb,&left_txs);
  av1_block_mode(aboveleft_mi, &aboveleft, &al_angle,&al_q,&al_sb,&al_txs);
  FILE *fp = fopen("intra_data.csv", "a");
  if (fp) {
    fprintf(fp, "%d,"
            "%d,%d,%d,"
            "%d,%d,%d,"
            "%d,%d,%d,"
            "%d,%d,%d,"
            "%d,%d,%d,",
            mode,
            above, left, aboveleft,
            above_angle, left_angle, al_angle,
            above_q, left_q, al_q,
            above_sb, left_sb, al_sb,
            above_txs, left_txs, al_txs
            );
    fprintf(fp, "\n");
    fclose(fp);
  }*/
}

static void write_inter_mode(aom_writer *w, PREDICTION_MODE mode,
                             FRAME_CONTEXT *ec_ctx, const int16_t mode_ctx) {
  const int16_t newmv_ctx = mode_ctx & NEWMV_CTX_MASK;

  aom_write_symbol(w, mode != NEWMV, ec_ctx->newmv_cdf[newmv_ctx], 2);

  if (mode != NEWMV) {
    const int16_t zeromv_ctx =
        (mode_ctx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
    aom_write_symbol(w, mode != GLOBALMV, ec_ctx->zeromv_cdf[zeromv_ctx], 2);

    if (mode != GLOBALMV) {
      int16_t refmv_ctx = (mode_ctx >> REFMV_OFFSET) & REFMV_CTX_MASK;
      aom_write_symbol(w, mode != NEARESTMV, ec_ctx->refmv_cdf[refmv_ctx], 2);
    }
  }
}

static void write_drl_idx(FRAME_CONTEXT *ec_ctx, const MB_MODE_INFO *mbmi,
                          const MB_MODE_INFO_EXT *mbmi_ext, aom_writer *w) {
  uint8_t ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);

  assert(mbmi->ref_mv_idx < 3);

  const int new_mv = mbmi->mode == NEWMV || mbmi->mode == NEW_NEWMV;
  if (new_mv) {
    int idx;
    for (idx = 0; idx < 2; ++idx) {
      if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx = av1_drl_ctx(mbmi_ext->weight[ref_frame_type], idx);

        aom_write_symbol(w, mbmi->ref_mv_idx != idx, ec_ctx->drl_cdf[drl_ctx],
                         2);
        if (mbmi->ref_mv_idx == idx) return;
      }
    }
    return;
  }

  if (have_nearmv_in_inter_mode(mbmi->mode)) {
    int idx;
    // TODO(jingning): Temporary solution to compensate the NEARESTMV offset.
    for (idx = 1; idx < 3; ++idx) {
      if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx = av1_drl_ctx(mbmi_ext->weight[ref_frame_type], idx);
        aom_write_symbol(w, mbmi->ref_mv_idx != (idx - 1),
                         ec_ctx->drl_cdf[drl_ctx], 2);
        if (mbmi->ref_mv_idx == (idx - 1)) return;
      }
    }
    return;
  }
}

static void write_inter_compound_mode(MACROBLOCKD *xd, aom_writer *w,
                                      PREDICTION_MODE mode,
                                      const int16_t mode_ctx) {
  assert(is_inter_compound_mode(mode));
  aom_write_symbol(w, INTER_COMPOUND_OFFSET(mode),
                   xd->tile_ctx->inter_compound_mode_cdf[mode_ctx],
                   INTER_COMPOUND_MODES);
}

static void write_tx_size_vartx(MACROBLOCKD *xd, const MB_MODE_INFO *mbmi,
                                TX_SIZE tx_size, int depth, int blk_row,
                                int blk_col, aom_writer *w) {
  FRAME_CONTEXT *const ec_ctx = xd->tile_ctx;
  const int max_blocks_high = max_block_high(xd, mbmi->sb_type, 0);
  const int max_blocks_wide = max_block_wide(xd, mbmi->sb_type, 0);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

  if (depth == MAX_VARTX_DEPTH) {
    txfm_partition_update(xd->above_txfm_context + blk_col,
                          xd->left_txfm_context + blk_row, tx_size, tx_size);
    return;
  }

  const int ctx = txfm_partition_context(xd->above_txfm_context + blk_col,
                                         xd->left_txfm_context + blk_row,
                                         mbmi->sb_type, tx_size);
  const int txb_size_index =
      av1_get_txb_size_index(mbmi->sb_type, blk_row, blk_col);
  const int write_txfm_partition =
      tx_size == mbmi->inter_tx_size[txb_size_index];
  if (write_txfm_partition) {
    aom_write_symbol(w, 0, ec_ctx->txfm_partition_cdf[ctx], 2);

    txfm_partition_update(xd->above_txfm_context + blk_col,
                          xd->left_txfm_context + blk_row, tx_size, tx_size);
    // TODO(yuec): set correct txfm partition update for qttx
  } else {
    const TX_SIZE sub_txs = sub_tx_size_map[tx_size];
    const int bsw = tx_size_wide_unit[sub_txs];
    const int bsh = tx_size_high_unit[sub_txs];

    aom_write_symbol(w, 1, ec_ctx->txfm_partition_cdf[ctx], 2);

    if (sub_txs == TX_4X4) {
      txfm_partition_update(xd->above_txfm_context + blk_col,
                            xd->left_txfm_context + blk_row, sub_txs, tx_size);
      return;
    }

    assert(bsw > 0 && bsh > 0);
    for (int row = 0; row < tx_size_high_unit[tx_size]; row += bsh)
      for (int col = 0; col < tx_size_wide_unit[tx_size]; col += bsw) {
        int offsetr = blk_row + row;
        int offsetc = blk_col + col;
        write_tx_size_vartx(xd, mbmi, sub_txs, depth + 1, offsetr, offsetc, w);
      }
  }
}

static void write_selected_tx_size(const MACROBLOCKD *xd, aom_writer *w) {
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  const BLOCK_SIZE bsize = mbmi->sb_type;
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;
  if (block_signals_txsize(bsize)) {
    const TX_SIZE tx_size = mbmi->tx_size;
    const int tx_size_ctx = get_tx_size_context(xd);
    const int depth = tx_size_to_depth(tx_size, bsize);
    const int max_depths = bsize_to_max_depth(bsize);
    const int32_t tx_size_cat = bsize_to_tx_size_cat(bsize);

    assert(depth >= 0 && depth <= max_depths);
    assert(!is_inter_block(mbmi));
    assert(IMPLIES(is_rect_tx(tx_size), is_rect_tx_allowed(xd, mbmi)));

    aom_write_symbol(w, depth, ec_ctx->tx_size_cdf[tx_size_cat][tx_size_ctx],
                     max_depths + 1);
  }
}

static int write_skip(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                      int segment_id, const MB_MODE_INFO *mi, aom_writer *w) {
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 1;
  } else {
    const int skip = mi->skip;
    const int ctx = av1_get_skip_context(xd);
    FRAME_CONTEXT *ec_ctx = xd->tile_ctx;
    aom_write_symbol(w, skip, ec_ctx->skip_cdfs[ctx], 2);
    return skip;
  }
}

static int write_skip_mode(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                           int segment_id, const MB_MODE_INFO *mi,
                           aom_writer *w) {
  if (!cm->current_frame.skip_mode_info.skip_mode_flag) return 0;
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP)) {
    return 0;
  }
  const int skip_mode = mi->skip_mode;
  if (!is_comp_ref_allowed(mi->sb_type)) {
    assert(!skip_mode);
    return 0;
  }
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME) ||
      segfeature_active(&cm->seg, segment_id, SEG_LVL_GLOBALMV)) {
    // These features imply single-reference mode, while skip mode implies
    // compound reference. Hence, the two are mutually exclusive.
    // In other words, skip_mode is implicitly 0 here.
    assert(!skip_mode);
    return 0;
  }
  const int ctx = av1_get_skip_mode_context(xd);
  aom_write_symbol(w, skip_mode, xd->tile_ctx->skip_mode_cdfs[ctx], 2);
  return skip_mode;
}

static void write_is_inter(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                           int segment_id, aom_writer *w, const int is_inter) {
  if (!segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    if (segfeature_active(&cm->seg, segment_id, SEG_LVL_GLOBALMV)) {
      assert(is_inter);
      return;
    }
    const int ctx = av1_get_intra_inter_context(xd);
    FRAME_CONTEXT *ec_ctx = xd->tile_ctx;
    aom_write_symbol(w, is_inter, ec_ctx->intra_inter_cdf[ctx], 2);
  }
}

static void write_motion_mode(const AV1_COMMON *cm, MACROBLOCKD *xd,
                              const MB_MODE_INFO *mbmi, aom_writer *w) {
  MOTION_MODE last_motion_mode_allowed =
      cm->switchable_motion_mode
          ? motion_mode_allowed(cm->global_motion, xd, mbmi,
                                cm->allow_warped_motion)
          : SIMPLE_TRANSLATION;
  assert(mbmi->motion_mode <= last_motion_mode_allowed);
  switch (last_motion_mode_allowed) {
    case SIMPLE_TRANSLATION: break;
    case OBMC_CAUSAL:
      aom_write_symbol(w, mbmi->motion_mode == OBMC_CAUSAL,
                       xd->tile_ctx->obmc_cdf[mbmi->sb_type], 2);
      break;
    default:
      aom_write_symbol(w, mbmi->motion_mode,
                       xd->tile_ctx->motion_mode_cdf[mbmi->sb_type],
                       MOTION_MODES);
  }
}

static void write_delta_qindex(const MACROBLOCKD *xd, int delta_qindex,
                               aom_writer *w) {
  int sign = delta_qindex < 0;
  int abs = sign ? -delta_qindex : delta_qindex;
  int rem_bits, thr;
  int smallval = abs < DELTA_Q_SMALL ? 1 : 0;
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;

  aom_write_symbol(w, AOMMIN(abs, DELTA_Q_SMALL), ec_ctx->delta_q_cdf,
                   DELTA_Q_PROBS + 1);

  if (!smallval) {
    rem_bits = get_msb(abs - 1);
    thr = (1 << rem_bits) + 1;
    aom_write_literal(w, rem_bits - 1, 3);
    aom_write_literal(w, abs - thr, rem_bits);
  }
  if (abs > 0) {
    aom_write_bit(w, sign);
  }
}

static void write_delta_lflevel(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                                int lf_id, int delta_lflevel, aom_writer *w) {
  int sign = delta_lflevel < 0;
  int abs = sign ? -delta_lflevel : delta_lflevel;
  int rem_bits, thr;
  int smallval = abs < DELTA_LF_SMALL ? 1 : 0;
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;

  if (cm->delta_q_info.delta_lf_multi) {
    assert(lf_id >= 0 && lf_id < (av1_num_planes(cm) > 1 ? FRAME_LF_COUNT
                                                         : FRAME_LF_COUNT - 2));
    aom_write_symbol(w, AOMMIN(abs, DELTA_LF_SMALL),
                     ec_ctx->delta_lf_multi_cdf[lf_id], DELTA_LF_PROBS + 1);
  } else {
    aom_write_symbol(w, AOMMIN(abs, DELTA_LF_SMALL), ec_ctx->delta_lf_cdf,
                     DELTA_LF_PROBS + 1);
  }

  if (!smallval) {
    rem_bits = get_msb(abs - 1);
    thr = (1 << rem_bits) + 1;
    aom_write_literal(w, rem_bits - 1, 3);
    aom_write_literal(w, abs - thr, rem_bits);
  }
  if (abs > 0) {
    aom_write_bit(w, sign);
  }
}

static void pack_map_tokens(aom_writer *w, const TOKENEXTRA **tp, int n,
                            int num) {
  const TOKENEXTRA *p = *tp;
  write_uniform(w, n, p->token);  // The first color index.
  ++p;
  --num;
  for (int i = 0; i < num; ++i) {
    aom_write_symbol(w, p->token, p->color_map_cdf, n);
    ++p;
  }
  *tp = p;
}

static void pack_txb_tokens(aom_writer *w, AV1_COMMON *cm, MACROBLOCK *const x,
                            const TOKENEXTRA **tp,
                            const TOKENEXTRA *const tok_end, MACROBLOCKD *xd,
                            MB_MODE_INFO *mbmi, int plane,
                            BLOCK_SIZE plane_bsize, aom_bit_depth_t bit_depth,
                            int block, int blk_row, int blk_col,
                            TX_SIZE tx_size, TOKEN_STATS *token_stats) {
  const int max_blocks_high = max_block_high(xd, plane_bsize, plane);
  const int max_blocks_wide = max_block_wide(xd, plane_bsize, plane);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const TX_SIZE plane_tx_size =
      plane ? av1_get_max_uv_txsize(mbmi->sb_type, pd->subsampling_x,
                                    pd->subsampling_y)
            : mbmi->inter_tx_size[av1_get_txb_size_index(plane_bsize, blk_row,
                                                         blk_col)];

  if (tx_size == plane_tx_size || plane) {
    const CB_COEFF_BUFFER *cb_coef_buff = x->cb_coef_buff;
    const int txb_offset =
        x->mbmi_ext->cb_offset / (TX_SIZE_W_MIN * TX_SIZE_H_MIN);
    const tran_low_t *tcoeff_txb =
        cb_coef_buff->tcoeff[plane] + x->mbmi_ext->cb_offset;
    const uint16_t *eob_txb = cb_coef_buff->eobs[plane] + txb_offset;
    const uint8_t *txb_skip_ctx_txb =
        cb_coef_buff->txb_skip_ctx[plane] + txb_offset;
    const int *dc_sign_ctx_txb = cb_coef_buff->dc_sign_ctx[plane] + txb_offset;
    const tran_low_t *tcoeff = BLOCK_OFFSET(tcoeff_txb, block);
    const uint16_t eob = eob_txb[block];
    TXB_CTX txb_ctx = { txb_skip_ctx_txb[block], dc_sign_ctx_txb[block] };
    av1_write_coeffs_txb(cm, xd, w, blk_row, blk_col, plane, tx_size, tcoeff,
                         eob, &txb_ctx);
#if CONFIG_RD_DEBUG
    TOKEN_STATS tmp_token_stats;
    init_token_stats(&tmp_token_stats);
    token_stats->txb_coeff_cost_map[blk_row][blk_col] = tmp_token_stats.cost;
    token_stats->cost += tmp_token_stats.cost;
#endif
  } else {
    const TX_SIZE sub_txs = sub_tx_size_map[tx_size];
    const int bsw = tx_size_wide_unit[sub_txs];
    const int bsh = tx_size_high_unit[sub_txs];
    const int step = bsh * bsw;

    assert(bsw > 0 && bsh > 0);

    for (int r = 0; r < tx_size_high_unit[tx_size]; r += bsh) {
      for (int c = 0; c < tx_size_wide_unit[tx_size]; c += bsw) {
        const int offsetr = blk_row + r;
        const int offsetc = blk_col + c;
        if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide) continue;
        pack_txb_tokens(w, cm, x, tp, tok_end, xd, mbmi, plane, plane_bsize,
                        bit_depth, block, offsetr, offsetc, sub_txs,
                        token_stats);
        block += step;
      }
    }
  }
}

static INLINE void set_spatial_segment_id(const AV1_COMMON *const cm,
                                          uint8_t *segment_ids,
                                          BLOCK_SIZE bsize, int mi_row,
                                          int mi_col, int segment_id) {
  const int mi_offset = mi_row * cm->mi_cols + mi_col;
  const int bw = mi_size_wide[bsize];
  const int bh = mi_size_high[bsize];
  const int xmis = AOMMIN(cm->mi_cols - mi_col, bw);
  const int ymis = AOMMIN(cm->mi_rows - mi_row, bh);
  int x, y;

  for (y = 0; y < ymis; ++y)
    for (x = 0; x < xmis; ++x)
      segment_ids[mi_offset + y * cm->mi_cols + x] = segment_id;
}

int av1_neg_interleave(int x, int ref, int max) {
  assert(x < max);
  const int diff = x - ref;
  if (!ref) return x;
  if (ref >= (max - 1)) return -x + max - 1;
  if (2 * ref < max) {
    if (abs(diff) <= ref) {
      if (diff > 0)
        return (diff << 1) - 1;
      else
        return ((-diff) << 1);
    }
    return x;
  } else {
    if (abs(diff) < (max - ref)) {
      if (diff > 0)
        return (diff << 1) - 1;
      else
        return ((-diff) << 1);
    }
    return (max - x) - 1;
  }
}

static void write_segment_id(AV1_COMP *cpi, const MB_MODE_INFO *const mbmi,
                             aom_writer *w, const struct segmentation *seg,
                             struct segmentation_probs *segp, int mi_row,
                             int mi_col, int skip) {
  if (!seg->enabled || !seg->update_map) return;

  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  int cdf_num;
  const int pred = av1_get_spatial_seg_pred(cm, xd, mi_row, mi_col, &cdf_num);

  if (skip) {
    // Still need to transmit tx size for intra blocks even if skip is
    // true. Changing segment_id may make the tx size become invalid, e.g
    // changing from lossless to lossy.
    assert(is_inter_block(mbmi) || !cpi->has_lossless_segment);

    set_spatial_segment_id(cm, cm->cur_frame->seg_map, mbmi->sb_type, mi_row,
                           mi_col, pred);
    set_spatial_segment_id(cm, cpi->segmentation_map, mbmi->sb_type, mi_row,
                           mi_col, pred);
    /* mbmi is read only but we need to update segment_id */
    ((MB_MODE_INFO *)mbmi)->segment_id = pred;
    return;
  }

  const int coded_id =
      av1_neg_interleave(mbmi->segment_id, pred, seg->last_active_segid + 1);
  aom_cdf_prob *pred_cdf = segp->spatial_pred_seg_cdf[cdf_num];
  aom_write_symbol(w, coded_id, pred_cdf, MAX_SEGMENTS);
  set_spatial_segment_id(cm, cm->cur_frame->seg_map, mbmi->sb_type, mi_row,
                         mi_col, mbmi->segment_id);
}

#define WRITE_REF_BIT(bname, pname) \
  aom_write_symbol(w, bname, av1_get_pred_cdf_##pname(xd), 2)

// This function encodes the reference frame
static void write_ref_frames(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                             aom_writer *w) {
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  const int is_compound = has_second_ref(mbmi);
  const int segment_id = mbmi->segment_id;

  // If segment level coding of this signal is disabled...
  // or the segment allows multiple reference frame options
  if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] ==
           get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME));
  } else if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP) ||
             segfeature_active(&cm->seg, segment_id, SEG_LVL_GLOBALMV)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] == LAST_FRAME);
  } else {
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
    if (cm->current_frame.reference_mode == REFERENCE_MODE_SELECT) {
      if (is_comp_ref_allowed(mbmi->sb_type))
        aom_write_symbol(w, is_compound, av1_get_reference_mode_cdf(xd), 2);
    } else {
      assert((!is_compound) ==
             (cm->current_frame.reference_mode == SINGLE_REFERENCE));
    }

    if (is_compound) {
      const COMP_REFERENCE_TYPE comp_ref_type = has_uni_comp_refs(mbmi)
                                                    ? UNIDIR_COMP_REFERENCE
                                                    : BIDIR_COMP_REFERENCE;
      aom_write_symbol(w, comp_ref_type, av1_get_comp_reference_type_cdf(xd),
                       2);

      if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
        const int bit = mbmi->ref_frame[0] == BWDREF_FRAME;
        WRITE_REF_BIT(bit, uni_comp_ref_p);

        if (!bit) {
          assert(mbmi->ref_frame[0] == LAST_FRAME);
          const int bit1 = mbmi->ref_frame[1] == LAST3_FRAME ||
                           mbmi->ref_frame[1] == GOLDEN_FRAME;
          WRITE_REF_BIT(bit1, uni_comp_ref_p1);
          if (bit1) {
            const int bit2 = mbmi->ref_frame[1] == GOLDEN_FRAME;
            WRITE_REF_BIT(bit2, uni_comp_ref_p2);
          }
        } else {
          assert(mbmi->ref_frame[1] == ALTREF_FRAME);
        }

        return;
      }

      assert(comp_ref_type == BIDIR_COMP_REFERENCE);

      const int bit = (mbmi->ref_frame[0] == GOLDEN_FRAME ||
                       mbmi->ref_frame[0] == LAST3_FRAME);
      WRITE_REF_BIT(bit, comp_ref_p);

      if (!bit) {
        const int bit1 = mbmi->ref_frame[0] == LAST2_FRAME;
        WRITE_REF_BIT(bit1, comp_ref_p1);
      } else {
        const int bit2 = mbmi->ref_frame[0] == GOLDEN_FRAME;
        WRITE_REF_BIT(bit2, comp_ref_p2);
      }

      const int bit_bwd = mbmi->ref_frame[1] == ALTREF_FRAME;
      WRITE_REF_BIT(bit_bwd, comp_bwdref_p);

      if (!bit_bwd) {
        WRITE_REF_BIT(mbmi->ref_frame[1] == ALTREF2_FRAME, comp_bwdref_p1);
      }

    } else {
      const int bit0 = (mbmi->ref_frame[0] <= ALTREF_FRAME &&
                        mbmi->ref_frame[0] >= BWDREF_FRAME);
      WRITE_REF_BIT(bit0, single_ref_p1);

      if (bit0) {
        const int bit1 = mbmi->ref_frame[0] == ALTREF_FRAME;
        WRITE_REF_BIT(bit1, single_ref_p2);

        if (!bit1) {
          WRITE_REF_BIT(mbmi->ref_frame[0] == ALTREF2_FRAME, single_ref_p6);
        }
      } else {
        const int bit2 = (mbmi->ref_frame[0] == LAST3_FRAME ||
                          mbmi->ref_frame[0] == GOLDEN_FRAME);
        WRITE_REF_BIT(bit2, single_ref_p3);

        if (!bit2) {
          const int bit3 = mbmi->ref_frame[0] != LAST_FRAME;
          WRITE_REF_BIT(bit3, single_ref_p4);
        } else {
          const int bit4 = mbmi->ref_frame[0] != LAST3_FRAME;
          WRITE_REF_BIT(bit4, single_ref_p5);
        }
      }
    }
  }
}

static void write_filter_intra_mode_info(const AV1_COMMON *cm,
                                         const MACROBLOCKD *xd,
                                         const MB_MODE_INFO *const mbmi,
                                         aom_writer *w) {
  if (av1_filter_intra_allowed(cm, mbmi)) {
    aom_write_symbol(w, mbmi->filter_intra_mode_info.use_filter_intra,
                     xd->tile_ctx->filter_intra_cdfs[mbmi->sb_type], 2);
    if (mbmi->filter_intra_mode_info.use_filter_intra) {
      const FILTER_INTRA_MODE mode =
          mbmi->filter_intra_mode_info.filter_intra_mode;
      aom_write_symbol(w, mode, xd->tile_ctx->filter_intra_mode_cdf,
                       FILTER_INTRA_MODES);
    }
  }
}

static void write_angle_delta(aom_writer *w, int angle_delta,
                              aom_cdf_prob *cdf) {
  aom_write_symbol(w, angle_delta + MAX_ANGLE_DELTA, cdf,
                   2 * MAX_ANGLE_DELTA + 1);
}

static void write_mb_interp_filter(AV1_COMP *cpi, const MACROBLOCKD *xd,
                                   aom_writer *w) {
  AV1_COMMON *const cm = &cpi->common;
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;

  if (!av1_is_interp_needed(xd)) {
    int_interpfilters filters =
        av1_broadcast_interp_filter(av1_unswitchable_filter(cm->interp_filter));
    assert(mbmi->interp_filters.as_int == filters.as_int);
    (void)filters;
    return;
  }
  if (cm->interp_filter == SWITCHABLE) {
    int dir;
    for (dir = 0; dir < 2; ++dir) {
      const int ctx = av1_get_pred_context_switchable_interp(xd, dir);
      InterpFilter filter =
          av1_extract_interp_filter(mbmi->interp_filters, dir);
      aom_write_symbol(w, filter, ec_ctx->switchable_interp_cdf[ctx],
                       SWITCHABLE_FILTERS);
      ++cm->cur_frame->interp_filter_selected[filter];
      if (cm->seq_params.enable_dual_filter == 0) return;
    }
  }
}

// Transmit color values with delta encoding. Write the first value as
// literal, and the deltas between each value and the previous one. "min_val" is
// the smallest possible value of the deltas.
static void delta_encode_palette_colors(const int *colors, int num,
                                        int bit_depth, int min_val,
                                        aom_writer *w) {
  if (num <= 0) return;
  assert(colors[0] < (1 << bit_depth));
  aom_write_literal(w, colors[0], bit_depth);
  if (num == 1) return;
  int max_delta = 0;
  int deltas[PALETTE_MAX_SIZE];
  memset(deltas, 0, sizeof(deltas));
  for (int i = 1; i < num; ++i) {
    assert(colors[i] < (1 << bit_depth));
    const int delta = colors[i] - colors[i - 1];
    deltas[i - 1] = delta;
    assert(delta >= min_val);
    if (delta > max_delta) max_delta = delta;
  }
  const int min_bits = bit_depth - 3;
  int bits = AOMMAX(av1_ceil_log2(max_delta + 1 - min_val), min_bits);
  assert(bits <= bit_depth);
  int range = (1 << bit_depth) - colors[0] - min_val;
  aom_write_literal(w, bits - min_bits, 2);
  for (int i = 0; i < num - 1; ++i) {
    aom_write_literal(w, deltas[i] - min_val, bits);
    range -= deltas[i];
    bits = AOMMIN(bits, av1_ceil_log2(range));
  }
}

// Transmit luma palette color values. First signal if each color in the color
// cache is used. Those colors that are not in the cache are transmitted with
// delta encoding.
static void write_palette_colors_y(const MACROBLOCKD *const xd,
                                   const PALETTE_MODE_INFO *const pmi,
                                   int bit_depth, aom_writer *w) {
  const int n = pmi->palette_size[0];
  uint16_t color_cache[2 * PALETTE_MAX_SIZE];
  const int n_cache = av1_get_palette_cache(xd, 0, color_cache);
  int out_cache_colors[PALETTE_MAX_SIZE];
  uint8_t cache_color_found[2 * PALETTE_MAX_SIZE];
  const int n_out_cache =
      av1_index_color_cache(color_cache, n_cache, pmi->palette_colors, n,
                            cache_color_found, out_cache_colors);
  int n_in_cache = 0;
  for (int i = 0; i < n_cache && n_in_cache < n; ++i) {
    const int found = cache_color_found[i];
    aom_write_bit(w, found);
    n_in_cache += found;
  }
  assert(n_in_cache + n_out_cache == n);
  delta_encode_palette_colors(out_cache_colors, n_out_cache, bit_depth, 1, w);
}

// Write chroma palette color values. U channel is handled similarly to the luma
// channel. For v channel, either use delta encoding or transmit raw values
// directly, whichever costs less.
static void write_palette_colors_uv(const MACROBLOCKD *const xd,
                                    const PALETTE_MODE_INFO *const pmi,
                                    int bit_depth, aom_writer *w) {
  const int n = pmi->palette_size[1];
  const uint16_t *colors_u = pmi->palette_colors + PALETTE_MAX_SIZE;
  const uint16_t *colors_v = pmi->palette_colors + 2 * PALETTE_MAX_SIZE;
  // U channel colors.
  uint16_t color_cache[2 * PALETTE_MAX_SIZE];
  const int n_cache = av1_get_palette_cache(xd, 1, color_cache);
  int out_cache_colors[PALETTE_MAX_SIZE];
  uint8_t cache_color_found[2 * PALETTE_MAX_SIZE];
  const int n_out_cache = av1_index_color_cache(
      color_cache, n_cache, colors_u, n, cache_color_found, out_cache_colors);
  int n_in_cache = 0;
  for (int i = 0; i < n_cache && n_in_cache < n; ++i) {
    const int found = cache_color_found[i];
    aom_write_bit(w, found);
    n_in_cache += found;
  }
  delta_encode_palette_colors(out_cache_colors, n_out_cache, bit_depth, 0, w);

  // V channel colors. Don't use color cache as the colors are not sorted.
  const int max_val = 1 << bit_depth;
  int zero_count = 0, min_bits_v = 0;
  int bits_v =
      av1_get_palette_delta_bits_v(pmi, bit_depth, &zero_count, &min_bits_v);
  const int rate_using_delta =
      2 + bit_depth + (bits_v + 1) * (n - 1) - zero_count;
  const int rate_using_raw = bit_depth * n;
  if (rate_using_delta < rate_using_raw) {  // delta encoding
    assert(colors_v[0] < (1 << bit_depth));
    aom_write_bit(w, 1);
    aom_write_literal(w, bits_v - min_bits_v, 2);
    aom_write_literal(w, colors_v[0], bit_depth);
    for (int i = 1; i < n; ++i) {
      assert(colors_v[i] < (1 << bit_depth));
      if (colors_v[i] == colors_v[i - 1]) {  // No need to signal sign bit.
        aom_write_literal(w, 0, bits_v);
        continue;
      }
      const int delta = abs((int)colors_v[i] - colors_v[i - 1]);
      const int sign_bit = colors_v[i] < colors_v[i - 1];
      if (delta <= max_val - delta) {
        aom_write_literal(w, delta, bits_v);
        aom_write_bit(w, sign_bit);
      } else {
        aom_write_literal(w, max_val - delta, bits_v);
        aom_write_bit(w, !sign_bit);
      }
    }
  } else {  // Transmit raw values.
    aom_write_bit(w, 0);
    for (int i = 0; i < n; ++i) {
      assert(colors_v[i] < (1 << bit_depth));
      aom_write_literal(w, colors_v[i], bit_depth);
    }
  }
}

static void write_palette_mode_info(const AV1_COMMON *cm, const MACROBLOCKD *xd,
                                    const MB_MODE_INFO *const mbmi, int mi_row,
                                    int mi_col, aom_writer *w) {
  const int num_planes = av1_num_planes(cm);
  const BLOCK_SIZE bsize = mbmi->sb_type;
  assert(av1_allow_palette(cm->allow_screen_content_tools, bsize));
  const PALETTE_MODE_INFO *const pmi = &mbmi->palette_mode_info;
  const int bsize_ctx = av1_get_palette_bsize_ctx(bsize);

  if (mbmi->mode == DC_PRED) {
    const int n = pmi->palette_size[0];
    const int palette_y_mode_ctx = av1_get_palette_mode_ctx(xd);
    aom_write_symbol(
        w, n > 0,
        xd->tile_ctx->palette_y_mode_cdf[bsize_ctx][palette_y_mode_ctx], 2);
    if (n > 0) {
      aom_write_symbol(w, n - PALETTE_MIN_SIZE,
                       xd->tile_ctx->palette_y_size_cdf[bsize_ctx],
                       PALETTE_SIZES);
      write_palette_colors_y(xd, pmi, cm->seq_params.bit_depth, w);
    }
  }

  const int uv_dc_pred =
      num_planes > 1 && mbmi->uv_mode == UV_DC_PRED &&
      is_chroma_reference(mi_row, mi_col, bsize, xd->plane[1].subsampling_x,
                          xd->plane[1].subsampling_y);
  if (uv_dc_pred) {
    const int n = pmi->palette_size[1];
    const int palette_uv_mode_ctx = (pmi->palette_size[0] > 0);
    aom_write_symbol(w, n > 0,
                     xd->tile_ctx->palette_uv_mode_cdf[palette_uv_mode_ctx], 2);
    if (n > 0) {
      aom_write_symbol(w, n - PALETTE_MIN_SIZE,
                       xd->tile_ctx->palette_uv_size_cdf[bsize_ctx],
                       PALETTE_SIZES);
      write_palette_colors_uv(xd, pmi, cm->seq_params.bit_depth, w);
    }
  }
}

void av1_write_tx_type(const AV1_COMMON *const cm, const MACROBLOCKD *xd,
                       int blk_row, int blk_col, int plane, TX_SIZE tx_size,
                       aom_writer *w) {
  MB_MODE_INFO *mbmi = xd->mi[0];
  const int is_inter = is_inter_block(mbmi);
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;

  // Only y plane's tx_type is transmitted
  if (plane > 0) return;
  PLANE_TYPE plane_type = get_plane_type(plane);
  TX_TYPE tx_type = av1_get_tx_type(plane_type, xd, blk_row, blk_col, tx_size,
                                    cm->reduced_tx_set_used);

  const TX_SIZE square_tx_size = txsize_sqr_map[tx_size];
  if (get_ext_tx_types(tx_size, is_inter, cm->reduced_tx_set_used) > 1 &&
      ((!cm->seg.enabled && cm->base_qindex > 0) ||
       (cm->seg.enabled && xd->qindex[mbmi->segment_id] > 0)) &&
      !mbmi->skip &&
      !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
    const TxSetType tx_set_type =
        av1_get_ext_tx_set_type(tx_size, is_inter, cm->reduced_tx_set_used);
    const int eset = get_ext_tx_set(tx_size, is_inter, cm->reduced_tx_set_used);
    // eset == 0 should correspond to a set with only DCT_DCT and there
    // is no need to send the tx_type
    assert(eset > 0);
    assert(av1_ext_tx_used[tx_set_type][tx_type]);
    if (is_inter) {
      aom_write_symbol(w, av1_ext_tx_ind[tx_set_type][tx_type],
                       ec_ctx->inter_ext_tx_cdf[eset][square_tx_size],
                       av1_num_ext_tx_set[tx_set_type]);
    } else {
      PREDICTION_MODE intra_dir;
      if (mbmi->filter_intra_mode_info.use_filter_intra)
        intra_dir =
            fimode_to_intradir[mbmi->filter_intra_mode_info.filter_intra_mode];
      else
        intra_dir = mbmi->mode;
      aom_write_symbol(
          w, av1_ext_tx_ind[tx_set_type][tx_type],
          ec_ctx->intra_ext_tx_cdf[eset][square_tx_size][intra_dir],
          av1_num_ext_tx_set[tx_set_type]);
    }
  }
}

static void write_intra_y_mode_nonkf(FRAME_CONTEXT *frame_ctx, BLOCK_SIZE bsize,
                                     PREDICTION_MODE mode, aom_writer *w) {
  aom_write_symbol(w, mode, frame_ctx->y_mode_cdf[size_group_lookup[bsize]],
                   INTRA_MODES);
}

static void write_intra_uv_mode(FRAME_CONTEXT *frame_ctx,
                                UV_PREDICTION_MODE uv_mode,
                                PREDICTION_MODE y_mode,
                                CFL_ALLOWED_TYPE cfl_allowed, aom_writer *w) {
  aom_write_symbol(w, uv_mode, frame_ctx->uv_mode_cdf[cfl_allowed][y_mode],
                   UV_INTRA_MODES - !cfl_allowed);
}

static void write_cfl_alphas(FRAME_CONTEXT *const ec_ctx, uint8_t idx,
                             int8_t joint_sign, aom_writer *w) {
  aom_write_symbol(w, joint_sign, ec_ctx->cfl_sign_cdf, CFL_JOINT_SIGNS);
  // Magnitudes are only signaled for nonzero codes.
  if (CFL_SIGN_U(joint_sign) != CFL_SIGN_ZERO) {
    aom_cdf_prob *cdf_u = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_U(joint_sign)];
    aom_write_symbol(w, CFL_IDX_U(idx), cdf_u, CFL_ALPHABET_SIZE);
  }
  if (CFL_SIGN_V(joint_sign) != CFL_SIGN_ZERO) {
    aom_cdf_prob *cdf_v = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_V(joint_sign)];
    aom_write_symbol(w, CFL_IDX_V(idx), cdf_v, CFL_ALPHABET_SIZE);
  }
}

static void write_cdef(AV1_COMMON *cm, MACROBLOCKD *const xd, aom_writer *w,
                       int skip, int mi_col, int mi_row) {
  if (cm->coded_lossless || cm->allow_intrabc) return;

  const int m = ~((1 << (6 - MI_SIZE_LOG2)) - 1);
  const MB_MODE_INFO *mbmi =
      cm->mi_grid_visible[(mi_row & m) * cm->mi_stride + (mi_col & m)];
  // Initialise when at top left part of the superblock
  if (!(mi_row & (cm->seq_params.mib_size - 1)) &&
      !(mi_col & (cm->seq_params.mib_size - 1))) {  // Top left?
    xd->cdef_preset[0] = xd->cdef_preset[1] = xd->cdef_preset[2] =
        xd->cdef_preset[3] = -1;
  }

  // Emit CDEF param at first non-skip coding block
  const int mask = 1 << (6 - MI_SIZE_LOG2);
  const int index = cm->seq_params.sb_size == BLOCK_128X128
                        ? !!(mi_col & mask) + 2 * !!(mi_row & mask)
                        : 0;
  if (xd->cdef_preset[index] == -1 && !skip) {
    aom_write_literal(w, mbmi->cdef_strength, cm->cdef_info.cdef_bits);
    xd->cdef_preset[index] = mbmi->cdef_strength;
  }
}

static void write_inter_segment_id(AV1_COMP *cpi, aom_writer *w,
                                   const struct segmentation *const seg,
                                   struct segmentation_probs *const segp,
                                   int mi_row, int mi_col, int skip,
                                   int preskip) {
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  AV1_COMMON *const cm = &cpi->common;

  if (seg->update_map) {
    if (preskip) {
      if (!seg->segid_preskip) return;
    } else {
      if (seg->segid_preskip) return;
      if (skip) {
        write_segment_id(cpi, mbmi, w, seg, segp, mi_row, mi_col, 1);
        if (seg->temporal_update) mbmi->seg_id_predicted = 0;
        return;
      }
    }
    if (seg->temporal_update) {
      const int pred_flag = mbmi->seg_id_predicted;
      aom_cdf_prob *pred_cdf = av1_get_pred_cdf_seg_id(segp, xd);
      aom_write_symbol(w, pred_flag, pred_cdf, 2);
      if (!pred_flag) {
        write_segment_id(cpi, mbmi, w, seg, segp, mi_row, mi_col, 0);
      }
      if (pred_flag) {
        set_spatial_segment_id(cm, cm->cur_frame->seg_map, mbmi->sb_type,
                               mi_row, mi_col, mbmi->segment_id);
      }
    } else {
      write_segment_id(cpi, mbmi, w, seg, segp, mi_row, mi_col, 0);
    }
  }
}

// If delta q is present, writes delta_q index.
// Also writes delta_q loop filter levels, if present.
static void write_delta_q_params(AV1_COMP *cpi, const int mi_row,
                                 const int mi_col, int skip, aom_writer *w) {
  AV1_COMMON *const cm = &cpi->common;
  const DeltaQInfo *const delta_q_info = &cm->delta_q_info;

  if (delta_q_info->delta_q_present_flag) {
    MACROBLOCK *const x = &cpi->td.mb;
    MACROBLOCKD *const xd = &x->e_mbd;
    const MB_MODE_INFO *const mbmi = xd->mi[0];
    const BLOCK_SIZE bsize = mbmi->sb_type;
    const int super_block_upper_left =
        ((mi_row & (cm->seq_params.mib_size - 1)) == 0) &&
        ((mi_col & (cm->seq_params.mib_size - 1)) == 0);

    if ((bsize != cm->seq_params.sb_size || skip == 0) &&
        super_block_upper_left) {
      assert(mbmi->current_qindex > 0);
      const int reduced_delta_qindex =
          (mbmi->current_qindex - xd->current_qindex) /
          delta_q_info->delta_q_res;
      write_delta_qindex(xd, reduced_delta_qindex, w);
      xd->current_qindex = mbmi->current_qindex;
      if (delta_q_info->delta_lf_present_flag) {
        if (delta_q_info->delta_lf_multi) {
          const int frame_lf_count =
              av1_num_planes(cm) > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
          for (int lf_id = 0; lf_id < frame_lf_count; ++lf_id) {
            int reduced_delta_lflevel =
                (mbmi->delta_lf[lf_id] - xd->delta_lf[lf_id]) /
                delta_q_info->delta_lf_res;
            write_delta_lflevel(cm, xd, lf_id, reduced_delta_lflevel, w);
            xd->delta_lf[lf_id] = mbmi->delta_lf[lf_id];
          }
        } else {
          int reduced_delta_lflevel =
              (mbmi->delta_lf_from_base - xd->delta_lf_from_base) /
              delta_q_info->delta_lf_res;
          write_delta_lflevel(cm, xd, -1, reduced_delta_lflevel, w);
          xd->delta_lf_from_base = mbmi->delta_lf_from_base;
        }
      }
    }
  }
}

static void write_intra_prediction_modes(AV1_COMP *cpi, const int mi_row,
                                         const int mi_col, int is_keyframe,
                                         aom_writer *w) {
  const AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->td.mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  const PREDICTION_MODE mode = mbmi->mode;
  const BLOCK_SIZE bsize = mbmi->sb_type;

  // Y mode.
  if (is_keyframe) {
    const MB_MODE_INFO *const above_mi = xd->above_mbmi;
    const MB_MODE_INFO *const left_mi = xd->left_mbmi;
    const MB_MODE_INFO *const aboveleft_mi = xd->aboveleft_mbmi;  // MC 2019
    write_intra_y_mode_kf(ec_ctx, mbmi, above_mi, left_mi, aboveleft_mi, mode,
                          w);
  } else {
    write_intra_y_mode_nonkf(ec_ctx, bsize, mode, w);
  }

  // Y angle delta.
  const int use_angle_delta = av1_use_angle_delta(bsize);
  if (use_angle_delta && av1_is_directional_mode(mode)) {
    write_angle_delta(w, mbmi->angle_delta[PLANE_TYPE_Y],
                      ec_ctx->angle_delta_cdf[mode - V_PRED]);
  }

  // UV mode and UV angle delta.
  if (!cm->seq_params.monochrome &&
      is_chroma_reference(mi_row, mi_col, bsize, xd->plane[1].subsampling_x,
                          xd->plane[1].subsampling_y)) {
    const UV_PREDICTION_MODE uv_mode = mbmi->uv_mode;
    write_intra_uv_mode(ec_ctx, uv_mode, mode, is_cfl_allowed(xd), w);
    if (uv_mode == UV_CFL_PRED)
      write_cfl_alphas(ec_ctx, mbmi->cfl_alpha_idx, mbmi->cfl_alpha_signs, w);
    if (use_angle_delta && av1_is_directional_mode(get_uv_mode(uv_mode))) {
      write_angle_delta(w, mbmi->angle_delta[PLANE_TYPE_UV],
                        ec_ctx->angle_delta_cdf[uv_mode - V_PRED]);
    }
  }

  // Palette.
  if (av1_allow_palette(cm->allow_screen_content_tools, bsize)) {
    write_palette_mode_info(cm, xd, mbmi, mi_row, mi_col, w);
  }

  // Filter intra.
  write_filter_intra_mode_info(cm, xd, mbmi, w);
}

static void pack_inter_mode_mvs(AV1_COMP *cpi, const int mi_row,
                                const int mi_col, aom_writer *w) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->td.mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;
  const struct segmentation *const seg = &cm->seg;
  struct segmentation_probs *const segp = &ec_ctx->seg;
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  const MB_MODE_INFO_EXT *const mbmi_ext = x->mbmi_ext;
  const PREDICTION_MODE mode = mbmi->mode;
  const int segment_id = mbmi->segment_id;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int allow_hp = cm->allow_high_precision_mv;
  const int is_inter = is_inter_block(mbmi);
  const int is_compound = has_second_ref(mbmi);
  int ref;

  write_inter_segment_id(cpi, w, seg, segp, mi_row, mi_col, 0, 1);

  write_skip_mode(cm, xd, segment_id, mbmi, w);

  assert(IMPLIES(mbmi->skip_mode, mbmi->skip));
  const int skip =
      mbmi->skip_mode ? 1 : write_skip(cm, xd, segment_id, mbmi, w);

  write_inter_segment_id(cpi, w, seg, segp, mi_row, mi_col, skip, 0);

  write_cdef(cm, xd, w, skip, mi_col, mi_row);

  write_delta_q_params(cpi, mi_row, mi_col, skip, w);

  if (!mbmi->skip_mode) write_is_inter(cm, xd, mbmi->segment_id, w, is_inter);

  if (mbmi->skip_mode) return;

  if (!is_inter) {
    write_intra_prediction_modes(cpi, mi_row, mi_col, 0, w);
  } else {
    int16_t mode_ctx;

    av1_collect_neighbors_ref_counts(xd);

    write_ref_frames(cm, xd, w);

    mode_ctx =
        av1_mode_context_analyzer(mbmi_ext->mode_context, mbmi->ref_frame);

    // If segment skip is not enabled code the mode.
    if (!segfeature_active(seg, segment_id, SEG_LVL_SKIP)) {
      if (is_inter_compound_mode(mode))
        write_inter_compound_mode(xd, w, mode, mode_ctx);
      else if (is_inter_singleref_mode(mode))
        write_inter_mode(w, mode, ec_ctx, mode_ctx);

      if (mode == NEWMV || mode == NEW_NEWMV || have_nearmv_in_inter_mode(mode))
        write_drl_idx(ec_ctx, mbmi, mbmi_ext, w);
      else
        assert(mbmi->ref_mv_idx == 0);
    }

    if (mode == NEWMV || mode == NEW_NEWMV) {
      for (ref = 0; ref < 1 + is_compound; ++ref) {
        nmv_context *nmvc = &ec_ctx->nmvc;
        const int_mv ref_mv = av1_get_ref_mv(x, ref);
        av1_encode_mv(cpi, w, &mbmi->mv[ref].as_mv, &ref_mv.as_mv, nmvc,
                      allow_hp);
      }
    } else if (mode == NEAREST_NEWMV || mode == NEAR_NEWMV) {
      nmv_context *nmvc = &ec_ctx->nmvc;
      const int_mv ref_mv = av1_get_ref_mv(x, 1);
      av1_encode_mv(cpi, w, &mbmi->mv[1].as_mv, &ref_mv.as_mv, nmvc, allow_hp);
    } else if (mode == NEW_NEARESTMV || mode == NEW_NEARMV) {
      nmv_context *nmvc = &ec_ctx->nmvc;
      const int_mv ref_mv = av1_get_ref_mv(x, 0);
      av1_encode_mv(cpi, w, &mbmi->mv[0].as_mv, &ref_mv.as_mv, nmvc, allow_hp);
    }

    if (cpi->common.current_frame.reference_mode != COMPOUND_REFERENCE &&
        cpi->common.seq_params.enable_interintra_compound &&
        is_interintra_allowed(mbmi)) {
      const int interintra = mbmi->ref_frame[1] == INTRA_FRAME;
      const int bsize_group = size_group_lookup[bsize];
      aom_write_symbol(w, interintra, ec_ctx->interintra_cdf[bsize_group], 2);
      if (interintra) {
        aom_write_symbol(w, mbmi->interintra_mode,
                         ec_ctx->interintra_mode_cdf[bsize_group],
                         INTERINTRA_MODES);
        if (is_interintra_wedge_used(bsize)) {
          aom_write_symbol(w, mbmi->use_wedge_interintra,
                           ec_ctx->wedge_interintra_cdf[bsize], 2);
          if (mbmi->use_wedge_interintra) {
            aom_write_symbol(w, mbmi->interintra_wedge_index,
                             ec_ctx->wedge_idx_cdf[bsize], 16);
          }
        }
      }
    }

    if (mbmi->ref_frame[1] != INTRA_FRAME) write_motion_mode(cm, xd, mbmi, w);

    // First write idx to indicate current compound inter prediction mode group
    // Group A (0): dist_wtd_comp, compound_average
    // Group B (1): interintra, compound_diffwtd, wedge
    if (has_second_ref(mbmi)) {
      const int masked_compound_used = is_any_masked_compound_used(bsize) &&
                                       cm->seq_params.enable_masked_compound;

      if (masked_compound_used) {
        const int ctx_comp_group_idx = get_comp_group_idx_context(xd);
        aom_write_symbol(w, mbmi->comp_group_idx,
                         ec_ctx->comp_group_idx_cdf[ctx_comp_group_idx], 2);
      } else {
        assert(mbmi->comp_group_idx == 0);
      }

      if (mbmi->comp_group_idx == 0) {
        if (mbmi->compound_idx)
          assert(mbmi->interinter_comp.type == COMPOUND_AVERAGE);

        if (cm->seq_params.order_hint_info.enable_dist_wtd_comp) {
          const int comp_index_ctx = get_comp_index_context(cm, xd);
          aom_write_symbol(w, mbmi->compound_idx,
                           ec_ctx->compound_index_cdf[comp_index_ctx], 2);
        } else {
          assert(mbmi->compound_idx == 1);
        }
      } else {
        assert(cpi->common.current_frame.reference_mode != SINGLE_REFERENCE &&
               is_inter_compound_mode(mbmi->mode) &&
               mbmi->motion_mode == SIMPLE_TRANSLATION);
        assert(masked_compound_used);
        // compound_diffwtd, wedge
        assert(mbmi->interinter_comp.type == COMPOUND_WEDGE ||
               mbmi->interinter_comp.type == COMPOUND_DIFFWTD);

        if (is_interinter_compound_used(COMPOUND_WEDGE, bsize))
          aom_write_symbol(w, mbmi->interinter_comp.type - COMPOUND_WEDGE,
                           ec_ctx->compound_type_cdf[bsize],
                           MASKED_COMPOUND_TYPES);

        if (mbmi->interinter_comp.type == COMPOUND_WEDGE) {
          assert(is_interinter_compound_used(COMPOUND_WEDGE, bsize));
          aom_write_symbol(w, mbmi->interinter_comp.wedge_index,
                           ec_ctx->wedge_idx_cdf[bsize], 16);
          aom_write_bit(w, mbmi->interinter_comp.wedge_sign);
        } else {
          assert(mbmi->interinter_comp.type == COMPOUND_DIFFWTD);
          aom_write_literal(w, mbmi->interinter_comp.mask_type,
                            MAX_DIFFWTD_MASK_BITS);
        }
      }
    }
    write_mb_interp_filter(cpi, xd, w);
  }
}

static void write_intrabc_info(MACROBLOCKD *xd,
                               const MB_MODE_INFO_EXT *mbmi_ext,
                               aom_writer *w) {
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  int use_intrabc = is_intrabc_block(mbmi);
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;
  aom_write_symbol(w, use_intrabc, ec_ctx->intrabc_cdf, 2);
  if (use_intrabc) {
    assert(mbmi->mode == DC_PRED);
    assert(mbmi->uv_mode == UV_DC_PRED);
    assert(mbmi->motion_mode == SIMPLE_TRANSLATION);
    int_mv dv_ref = mbmi_ext->ref_mv_stack[INTRA_FRAME][0].this_mv;
    av1_encode_dv(w, &mbmi->mv[0].as_mv, &dv_ref.as_mv, &ec_ctx->ndvc);
  }
}

static void write_mb_modes_kf(AV1_COMP *cpi, MACROBLOCKD *xd,
                              const MB_MODE_INFO_EXT *mbmi_ext,
                              const int mi_row, const int mi_col,
                              aom_writer *w) {
  AV1_COMMON *const cm = &cpi->common;
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;
  const struct segmentation *const seg = &cm->seg;
  struct segmentation_probs *const segp = &ec_ctx->seg;
  const MB_MODE_INFO *const mbmi = xd->mi[0];

  if (seg->segid_preskip && seg->update_map)
    write_segment_id(cpi, mbmi, w, seg, segp, mi_row, mi_col, 0);

  const int skip = write_skip(cm, xd, mbmi->segment_id, mbmi, w);

  if (!seg->segid_preskip && seg->update_map)
    write_segment_id(cpi, mbmi, w, seg, segp, mi_row, mi_col, skip);

  write_cdef(cm, xd, w, skip, mi_col, mi_row);

  write_delta_q_params(cpi, mi_row, mi_col, skip, w);

  if (av1_allow_intrabc(cm)) {
    write_intrabc_info(xd, mbmi_ext, w);
    if (is_intrabc_block(mbmi)) return;
  }

  write_intra_prediction_modes(cpi, mi_row, mi_col, 1, w);
}

#if CONFIG_RD_DEBUG
static void dump_mode_info(MB_MODE_INFO *mi) {
  printf("\nmi->mi_row == %d\n", mi->mi_row);
  printf("&& mi->mi_col == %d\n", mi->mi_col);
  printf("&& mi->sb_type == %d\n", mi->sb_type);
  printf("&& mi->tx_size == %d\n", mi->tx_size);
  printf("&& mi->mode == %d\n", mi->mode);
}

static int rd_token_stats_mismatch(RD_STATS *rd_stats, TOKEN_STATS *token_stats,
                                   int plane) {
  if (rd_stats->txb_coeff_cost[plane] != token_stats->cost) {
    int r, c;
    printf("\nplane %d rd_stats->txb_coeff_cost %d token_stats->cost %d\n",
           plane, rd_stats->txb_coeff_cost[plane], token_stats->cost);
    printf("rd txb_coeff_cost_map\n");
    for (r = 0; r < TXB_COEFF_COST_MAP_SIZE; ++r) {
      for (c = 0; c < TXB_COEFF_COST_MAP_SIZE; ++c) {
        printf("%d ", rd_stats->txb_coeff_cost_map[plane][r][c]);
      }
      printf("\n");
    }

    printf("pack txb_coeff_cost_map\n");
    for (r = 0; r < TXB_COEFF_COST_MAP_SIZE; ++r) {
      for (c = 0; c < TXB_COEFF_COST_MAP_SIZE; ++c) {
        printf("%d ", token_stats->txb_coeff_cost_map[r][c]);
      }
      printf("\n");
    }
    return 1;
  }
  return 0;
}
#endif

#if ENC_MISMATCH_DEBUG
static void enc_dump_logs(AV1_COMP *cpi, int mi_row, int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const MB_MODE_INFO *const *mbmi =
      *(cm->mi_grid_visible + (mi_row * cm->mi_stride + mi_col));
  const MB_MODE_INFO_EXT *const *mbmi_ext =
      cpi->mbmi_ext_base + (mi_row * cm->mi_cols + mi_col);
  if (is_inter_block(mbmi)) {
#define FRAME_TO_CHECK 11
    if (cm->current_frame.frame_number == FRAME_TO_CHECK &&
        cm->show_frame == 1) {
      const BLOCK_SIZE bsize = mbmi->sb_type;

      int_mv mv[2] = { 0 };
      const int is_comp_ref = has_second_ref(mbmi);

      for (int ref = 0; ref < 1 + is_comp_ref; ++ref)
        mv[ref].as_mv = mbmi->mv[ref].as_mv;

      if (!is_comp_ref) {
        mv[1].as_int = 0;
      }

      const int16_t mode_ctx =
          is_comp_ref ? 0
                      : av1_mode_context_analyzer(mbmi_ext->mode_context,
                                                  mbmi->ref_frame);

      const int16_t newmv_ctx = mode_ctx & NEWMV_CTX_MASK;
      int16_t zeromv_ctx = -1;
      int16_t refmv_ctx = -1;

      if (mbmi->mode != NEWMV) {
        zeromv_ctx = (mode_ctx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
        if (mbmi->mode != GLOBALMV)
          refmv_ctx = (mode_ctx >> REFMV_OFFSET) & REFMV_CTX_MASK;
      }

      printf(
          "=== ENCODER ===: "
          "Frame=%d, (mi_row,mi_col)=(%d,%d), skip_mode=%d, mode=%d, bsize=%d, "
          "show_frame=%d, mv[0]=(%d,%d), mv[1]=(%d,%d), ref[0]=%d, "
          "ref[1]=%d, motion_mode=%d, mode_ctx=%d, "
          "newmv_ctx=%d, zeromv_ctx=%d, refmv_ctx=%d, tx_size=%d\n",
          cm->current_frame.frame_number, mi_row, mi_col, mbmi->skip_mode,
          mbmi->mode, bsize, cm->show_frame, mv[0].as_mv.row, mv[0].as_mv.col,
          mv[1].as_mv.row, mv[1].as_mv.col, mbmi->ref_frame[0],
          mbmi->ref_frame[1], mbmi->motion_mode, mode_ctx, newmv_ctx,
          zeromv_ctx, refmv_ctx, mbmi->tx_size);
    }
  }
}
#endif  // ENC_MISMATCH_DEBUG

static void write_mbmi_b(AV1_COMP *cpi, aom_writer *w, int mi_row, int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  MB_MODE_INFO *m = xd->mi[0];

  if (frame_is_intra_only(cm)) {
    write_mb_modes_kf(cpi, xd, cpi->td.mb.mbmi_ext, mi_row, mi_col, w);
  } else {
    // has_subpel_mv_component needs the ref frame buffers set up to look
    // up if they are scaled. has_subpel_mv_component is in turn needed by
    // write_switchable_interp_filter, which is called by pack_inter_mode_mvs.
    set_ref_ptrs(cm, xd, m->ref_frame[0], m->ref_frame[1]);

#if ENC_MISMATCH_DEBUG
    enc_dump_logs(cpi, mi_row, mi_col);
#endif  // ENC_MISMATCH_DEBUG

    pack_inter_mode_mvs(cpi, mi_row, mi_col, w);
  }
}

static void write_inter_txb_coeff(AV1_COMMON *const cm, MACROBLOCK *const x,
                                  MB_MODE_INFO *const mbmi, aom_writer *w,
                                  const TOKENEXTRA **tok,
                                  const TOKENEXTRA *const tok_end,
                                  TOKEN_STATS *token_stats, const int row,
                                  const int col, int *block, const int plane) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const BLOCK_SIZE bsize = mbmi->sb_type;
  assert(bsize < BLOCK_SIZES_ALL);
  const BLOCK_SIZE bsizec =
      scale_chroma_bsize(bsize, pd->subsampling_x, pd->subsampling_y);

  const BLOCK_SIZE plane_bsize =
      get_plane_block_size(bsizec, pd->subsampling_x, pd->subsampling_y);

  const TX_SIZE max_tx_size = get_vartx_max_txsize(xd, plane_bsize, plane);
  const int step =
      tx_size_wide_unit[max_tx_size] * tx_size_high_unit[max_tx_size];
  const int bkw = tx_size_wide_unit[max_tx_size];
  const int bkh = tx_size_high_unit[max_tx_size];

  const BLOCK_SIZE max_unit_bsize =
      get_plane_block_size(BLOCK_64X64, pd->subsampling_x, pd->subsampling_y);
  int mu_blocks_wide = block_size_wide[max_unit_bsize] >> tx_size_wide_log2[0];
  int mu_blocks_high = block_size_high[max_unit_bsize] >> tx_size_high_log2[0];

  int blk_row, blk_col;

  assert(plane_bsize < BLOCK_SIZES_ALL);
  const int num_4x4_w = block_size_wide[plane_bsize] >> tx_size_wide_log2[0];
  const int num_4x4_h = block_size_high[plane_bsize] >> tx_size_high_log2[0];

  const int unit_height =
      AOMMIN(mu_blocks_high + (row >> pd->subsampling_y), num_4x4_h);
  const int unit_width =
      AOMMIN(mu_blocks_wide + (col >> pd->subsampling_x), num_4x4_w);
  for (blk_row = row >> pd->subsampling_y; blk_row < unit_height;
       blk_row += bkh) {
    for (blk_col = col >> pd->subsampling_x; blk_col < unit_width;
         blk_col += bkw) {
      pack_txb_tokens(w, cm, x, tok, tok_end, xd, mbmi, plane, plane_bsize,
                      cm->seq_params.bit_depth, *block, blk_row, blk_col,
                      max_tx_size, token_stats);
      *block += step;
    }
  }
}

static void write_tokens_b(AV1_COMP *cpi, aom_writer *w, const TOKENEXTRA **tok,
                           const TOKENEXTRA *const tok_end, int mi_row,
                           int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->td.mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  const BLOCK_SIZE bsize = mbmi->sb_type;

  assert(!mbmi->skip);

  const int is_inter = is_inter_block(mbmi);
  if (!is_inter) {
    av1_write_coeffs_mb(cm, x, mi_row, mi_col, w, bsize);
  } else {
    int block[MAX_MB_PLANE] = { 0 };
    assert(bsize == get_plane_block_size(bsize, xd->plane[0].subsampling_x,
                                         xd->plane[0].subsampling_y));
    const int num_4x4_w = block_size_wide[bsize] >> tx_size_wide_log2[0];
    const int num_4x4_h = block_size_high[bsize] >> tx_size_high_log2[0];
    TOKEN_STATS token_stats;
    init_token_stats(&token_stats);

    const BLOCK_SIZE max_unit_bsize = BLOCK_64X64;
    assert(max_unit_bsize == get_plane_block_size(BLOCK_64X64,
                                                  xd->plane[0].subsampling_x,
                                                  xd->plane[0].subsampling_y));
    int mu_blocks_wide =
        block_size_wide[max_unit_bsize] >> tx_size_wide_log2[0];
    int mu_blocks_high =
        block_size_high[max_unit_bsize] >> tx_size_high_log2[0];

    mu_blocks_wide = AOMMIN(num_4x4_w, mu_blocks_wide);
    mu_blocks_high = AOMMIN(num_4x4_h, mu_blocks_high);

    const int num_planes = av1_num_planes(cm);
    for (int row = 0; row < num_4x4_h; row += mu_blocks_high) {
      for (int col = 0; col < num_4x4_w; col += mu_blocks_wide) {
        for (int plane = 0; plane < num_planes; ++plane) {
          const struct macroblockd_plane *const pd = &xd->plane[plane];
          if (!is_chroma_reference(mi_row, mi_col, bsize, pd->subsampling_x,
                                   pd->subsampling_y)) {
            continue;
          }
          write_inter_txb_coeff(cm, x, mbmi, w, tok, tok_end, &token_stats, row,
                                col, &block[plane], plane);
        }
      }
    }
#if CONFIG_RD_DEBUG
    for (int plane = 0; plane < num_planes; ++plane) {
      if (mbmi->sb_type >= BLOCK_8X8 &&
          rd_token_stats_mismatch(&mbmi->rd_stats, &token_stats, plane)) {
        dump_mode_info(mbmi);
        assert(0);
      }
    }
#endif  // CONFIG_RD_DEBUG
  }
}

static void write_modes_b(AV1_COMP *cpi, const TileInfo *const tile,
                          aom_writer *w, const TOKENEXTRA **tok,
                          const TOKENEXTRA *const tok_end, int mi_row,
                          int mi_col) {
  const AV1_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &cpi->td.mb.e_mbd;
  xd->mi = cm->mi_grid_visible + (mi_row * cm->mi_stride + mi_col);
  cpi->td.mb.mbmi_ext = cpi->mbmi_ext_base + (mi_row * cm->mi_cols + mi_col);

  const MB_MODE_INFO *mbmi = xd->mi[0];
  const BLOCK_SIZE bsize = mbmi->sb_type;
  assert(bsize <= cm->seq_params.sb_size ||
         (bsize >= BLOCK_SIZES && bsize < BLOCK_SIZES_ALL));

  const int bh = mi_size_high[bsize];
  const int bw = mi_size_wide[bsize];
  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);

  xd->above_txfm_context = cm->above_txfm_context[tile->tile_row] + mi_col;
  xd->left_txfm_context =
      xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);

  write_mbmi_b(cpi, w, mi_row, mi_col);

  for (int plane = 0; plane < AOMMIN(2, av1_num_planes(cm)); ++plane) {
    const uint8_t palette_size_plane =
        mbmi->palette_mode_info.palette_size[plane];
    assert(!mbmi->skip_mode || !palette_size_plane);
    if (palette_size_plane > 0) {
      assert(mbmi->use_intrabc == 0);
      assert(av1_allow_palette(cm->allow_screen_content_tools, mbmi->sb_type));
      int rows, cols;
      av1_get_block_dimensions(mbmi->sb_type, plane, xd, NULL, NULL, &rows,
                               &cols);
      assert(*tok < tok_end);
      pack_map_tokens(w, tok, palette_size_plane, rows * cols);
    }
  }

  const int is_inter_tx = is_inter_block(mbmi);
  const int skip = mbmi->skip;
  const int segment_id = mbmi->segment_id;
  if (cm->tx_mode == TX_MODE_SELECT && block_signals_txsize(bsize) &&
      !(is_inter_tx && skip) && !xd->lossless[segment_id]) {
    if (is_inter_tx) {  // This implies skip flag is 0.
      const TX_SIZE max_tx_size = get_vartx_max_txsize(xd, bsize, 0);
      const int txbh = tx_size_high_unit[max_tx_size];
      const int txbw = tx_size_wide_unit[max_tx_size];
      const int width = block_size_wide[bsize] >> tx_size_wide_log2[0];
      const int height = block_size_high[bsize] >> tx_size_high_log2[0];
      for (int idy = 0; idy < height; idy += txbh) {
        for (int idx = 0; idx < width; idx += txbw) {
          write_tx_size_vartx(xd, mbmi, max_tx_size, 0, idy, idx, w);
        }
      }
    } else {
      write_selected_tx_size(xd, w);
      set_txfm_ctxs(mbmi->tx_size, xd->n4_w, xd->n4_h, 0, xd);
    }
  } else {
    set_txfm_ctxs(mbmi->tx_size, xd->n4_w, xd->n4_h, skip && is_inter_tx, xd);
  }

  if (!mbmi->skip) {
    write_tokens_b(cpi, w, tok, tok_end, mi_row, mi_col);
  }
}

static void write_partition(const AV1_COMMON *const cm,
                            const MACROBLOCKD *const xd, int hbs, int mi_row,
                            int mi_col, PARTITION_TYPE p, BLOCK_SIZE bsize,
                            aom_writer *w) {
  const int is_partition_point = bsize >= BLOCK_8X8;

  if (!is_partition_point) return;

  const int has_rows = (mi_row + hbs) < cm->mi_rows;
  const int has_cols = (mi_col + hbs) < cm->mi_cols;
  const int ctx = partition_plane_context(xd, mi_row, mi_col, bsize);
  FRAME_CONTEXT *ec_ctx = xd->tile_ctx;

  if (!has_rows && !has_cols) {
    assert(p == PARTITION_SPLIT);
    return;
  }

  if (has_rows && has_cols) {
    aom_write_symbol(w, p, ec_ctx->partition_cdf[ctx],
                     partition_cdf_length(bsize));
  } else if (!has_rows && has_cols) {
    assert(p == PARTITION_SPLIT || p == PARTITION_HORZ);
    assert(bsize > BLOCK_8X8);
    aom_cdf_prob cdf[2];
    partition_gather_vert_alike(cdf, ec_ctx->partition_cdf[ctx], bsize);
    aom_write_cdf(w, p == PARTITION_SPLIT, cdf, 2);
  } else {
    assert(has_rows && !has_cols);
    assert(p == PARTITION_SPLIT || p == PARTITION_VERT);
    assert(bsize > BLOCK_8X8);
    aom_cdf_prob cdf[2];
    partition_gather_horz_alike(cdf, ec_ctx->partition_cdf[ctx], bsize);
    aom_write_cdf(w, p == PARTITION_SPLIT, cdf, 2);
  }
}

static void write_modes_sb(AV1_COMP *const cpi, const TileInfo *const tile,
                           aom_writer *const w, const TOKENEXTRA **tok,
                           const TOKENEXTRA *const tok_end, int mi_row,
                           int mi_col, BLOCK_SIZE bsize) {
  const AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  assert(bsize < BLOCK_SIZES_ALL);
  const int hbs = mi_size_wide[bsize] / 2;
  const int quarter_step = mi_size_wide[bsize] / 4;
  int i;
  const PARTITION_TYPE partition = get_partition(cm, mi_row, mi_col, bsize);
  const BLOCK_SIZE subsize = get_partition_subsize(bsize, partition);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  const int num_planes = av1_num_planes(cm);
  for (int plane = 0; plane < num_planes; ++plane) {
    int rcol0, rcol1, rrow0, rrow1;
    if (av1_loop_restoration_corners_in_sb(cm, plane, mi_row, mi_col, bsize,
                                           &rcol0, &rcol1, &rrow0, &rrow1)) {
      const int rstride = cm->rst_info[plane].horz_units_per_tile;
      for (int rrow = rrow0; rrow < rrow1; ++rrow) {
        for (int rcol = rcol0; rcol < rcol1; ++rcol) {
          const int runit_idx = rcol + rrow * rstride;
          const RestorationUnitInfo *rui =
              &cm->rst_info[plane].unit_info[runit_idx];
          loop_restoration_write_sb_coeffs(cm, xd, rui, w, plane,
                                           cpi->td.counts);
        }
      }
    }
  }

  write_partition(cm, xd, hbs, mi_row, mi_col, partition, bsize, w);
  switch (partition) {
    case PARTITION_NONE:
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col);
      break;
    case PARTITION_HORZ:
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col);
      if (mi_row + hbs < cm->mi_rows)
        write_modes_b(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col);
      break;
    case PARTITION_VERT:
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col);
      if (mi_col + hbs < cm->mi_cols)
        write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col + hbs);
      break;
    case PARTITION_SPLIT:
      write_modes_sb(cpi, tile, w, tok, tok_end, mi_row, mi_col, subsize);
      write_modes_sb(cpi, tile, w, tok, tok_end, mi_row, mi_col + hbs, subsize);
      write_modes_sb(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col, subsize);
      write_modes_sb(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col + hbs,
                     subsize);
      break;
    case PARTITION_HORZ_A:
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col + hbs);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col);
      break;
    case PARTITION_HORZ_B:
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col + hbs);
      break;
    case PARTITION_VERT_A:
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col + hbs);
      break;
    case PARTITION_VERT_B:
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row, mi_col + hbs);
      write_modes_b(cpi, tile, w, tok, tok_end, mi_row + hbs, mi_col + hbs);
      break;
    case PARTITION_HORZ_4:
      for (i = 0; i < 4; ++i) {
        int this_mi_row = mi_row + i * quarter_step;
        if (i > 0 && this_mi_row >= cm->mi_rows) break;

        write_modes_b(cpi, tile, w, tok, tok_end, this_mi_row, mi_col);
      }
      break;
    case PARTITION_VERT_4:
      for (i = 0; i < 4; ++i) {
        int this_mi_col = mi_col + i * quarter_step;
        if (i > 0 && this_mi_col >= cm->mi_cols) break;

        write_modes_b(cpi, tile, w, tok, tok_end, mi_row, this_mi_col);
      }
      break;
    default: assert(0);
  }

  // update partition context
  update_ext_partition_context(xd, mi_row, mi_col, subsize, bsize, partition);
}

static void write_modes(AV1_COMP *const cpi, const TileInfo *const tile,
                        aom_writer *const w, int tile_row, int tile_col) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  const int mi_row_start = tile->mi_row_start;
  const int mi_row_end = tile->mi_row_end;
  const int mi_col_start = tile->mi_col_start;
  const int mi_col_end = tile->mi_col_end;
  int mi_row, mi_col, sb_row_in_tile;

  av1_zero_above_context(cm, xd, mi_col_start, mi_col_end, tile->tile_row);
  av1_init_above_context(cm, xd, tile->tile_row);

  if (cpi->common.delta_q_info.delta_q_present_flag) {
    xd->current_qindex = cpi->common.base_qindex;
    if (cpi->common.delta_q_info.delta_lf_present_flag) {
      av1_reset_loop_filter_delta(xd, av1_num_planes(cm));
    }
  }

  for (mi_row = mi_row_start; mi_row < mi_row_end;
       mi_row += cm->seq_params.mib_size) {
    sb_row_in_tile =
        (mi_row - tile->mi_row_start) >> cm->seq_params.mib_size_log2;
    const TOKENEXTRA *tok =
        cpi->tplist[tile_row][tile_col][sb_row_in_tile].start;
    const TOKENEXTRA *tok_end =
        tok + cpi->tplist[tile_row][tile_col][sb_row_in_tile].count;

    av1_zero_left_context(xd);

    for (mi_col = mi_col_start; mi_col < mi_col_end;
         mi_col += cm->seq_params.mib_size) {
      cpi->td.mb.cb_coef_buff = av1_get_cb_coeff_buffer(cpi, mi_row, mi_col);
      write_modes_sb(cpi, tile, w, &tok, tok_end, mi_row, mi_col,
                     cm->seq_params.sb_size);
    }
    assert(tok == cpi->tplist[tile_row][tile_col][sb_row_in_tile].stop);
  }
}

static void encode_restoration_mode(AV1_COMMON *cm,
                                    struct aom_write_bit_buffer *wb) {
  assert(!cm->all_lossless);
  if (!cm->seq_params.enable_restoration) return;
  if (cm->allow_intrabc) return;
  const int num_planes = av1_num_planes(cm);
  int all_none = 1, chroma_none = 1;
  for (int p = 0; p < num_planes; ++p) {
    RestorationInfo *rsi = &cm->rst_info[p];
    if (rsi->frame_restoration_type != RESTORE_NONE) {
      all_none = 0;
      chroma_none &= p == 0;
    }
    switch (rsi->frame_restoration_type) {
      case RESTORE_NONE:
        aom_wb_write_bit(wb, 0);
        aom_wb_write_bit(wb, 0);
        break;
      case RESTORE_WIENER:
        aom_wb_write_bit(wb, 1);
        aom_wb_write_bit(wb, 0);
        break;
      case RESTORE_SGRPROJ:
        aom_wb_write_bit(wb, 1);
        aom_wb_write_bit(wb, 1);
        break;
      case RESTORE_SWITCHABLE:
        aom_wb_write_bit(wb, 0);
        aom_wb_write_bit(wb, 1);
        break;
      default: assert(0);
    }
  }
  if (!all_none) {
    assert(cm->seq_params.sb_size == BLOCK_64X64 ||
           cm->seq_params.sb_size == BLOCK_128X128);
    const int sb_size = cm->seq_params.sb_size == BLOCK_128X128 ? 128 : 64;

    RestorationInfo *rsi = &cm->rst_info[0];

    assert(rsi->restoration_unit_size >= sb_size);
    assert(RESTORATION_UNITSIZE_MAX == 256);

    if (sb_size == 64) {
      aom_wb_write_bit(wb, rsi->restoration_unit_size > 64);
    }
    if (rsi->restoration_unit_size > 64) {
      aom_wb_write_bit(wb, rsi->restoration_unit_size > 128);
    }
  }

  if (num_planes > 1) {
    int s = AOMMIN(cm->seq_params.subsampling_x, cm->seq_params.subsampling_y);
    if (s && !chroma_none) {
      aom_wb_write_bit(wb, cm->rst_info[1].restoration_unit_size !=
                               cm->rst_info[0].restoration_unit_size);
      assert(cm->rst_info[1].restoration_unit_size ==
                 cm->rst_info[0].restoration_unit_size ||
             cm->rst_info[1].restoration_unit_size ==
                 (cm->rst_info[0].restoration_unit_size >> s));
      assert(cm->rst_info[2].restoration_unit_size ==
             cm->rst_info[1].restoration_unit_size);
    } else if (!s) {
      assert(cm->rst_info[1].restoration_unit_size ==
             cm->rst_info[0].restoration_unit_size);
      assert(cm->rst_info[2].restoration_unit_size ==
             cm->rst_info[1].restoration_unit_size);
    }
  }
}

static void write_wiener_filter(int wiener_win, const WienerInfo *wiener_info,
                                WienerInfo *ref_wiener_info, aom_writer *wb) {
  if (wiener_win == WIENER_WIN)
    aom_write_primitive_refsubexpfin(
        wb, WIENER_FILT_TAP0_MAXV - WIENER_FILT_TAP0_MINV + 1,
        WIENER_FILT_TAP0_SUBEXP_K,
        ref_wiener_info->vfilter[0] - WIENER_FILT_TAP0_MINV,
        wiener_info->vfilter[0] - WIENER_FILT_TAP0_MINV);
  else
    assert(wiener_info->vfilter[0] == 0 &&
           wiener_info->vfilter[WIENER_WIN - 1] == 0);
  aom_write_primitive_refsubexpfin(
      wb, WIENER_FILT_TAP1_MAXV - WIENER_FILT_TAP1_MINV + 1,
      WIENER_FILT_TAP1_SUBEXP_K,
      ref_wiener_info->vfilter[1] - WIENER_FILT_TAP1_MINV,
      wiener_info->vfilter[1] - WIENER_FILT_TAP1_MINV);
  aom_write_primitive_refsubexpfin(
      wb, WIENER_FILT_TAP2_MAXV - WIENER_FILT_TAP2_MINV + 1,
      WIENER_FILT_TAP2_SUBEXP_K,
      ref_wiener_info->vfilter[2] - WIENER_FILT_TAP2_MINV,
      wiener_info->vfilter[2] - WIENER_FILT_TAP2_MINV);
  if (wiener_win == WIENER_WIN)
    aom_write_primitive_refsubexpfin(
        wb, WIENER_FILT_TAP0_MAXV - WIENER_FILT_TAP0_MINV + 1,
        WIENER_FILT_TAP0_SUBEXP_K,
        ref_wiener_info->hfilter[0] - WIENER_FILT_TAP0_MINV,
        wiener_info->hfilter[0] - WIENER_FILT_TAP0_MINV);
  else
    assert(wiener_info->hfilter[0] == 0 &&
           wiener_info->hfilter[WIENER_WIN - 1] == 0);
  aom_write_primitive_refsubexpfin(
      wb, WIENER_FILT_TAP1_MAXV - WIENER_FILT_TAP1_MINV + 1,
      WIENER_FILT_TAP1_SUBEXP_K,
      ref_wiener_info->hfilter[1] - WIENER_FILT_TAP1_MINV,
      wiener_info->hfilter[1] - WIENER_FILT_TAP1_MINV);
  aom_write_primitive_refsubexpfin(
      wb, WIENER_FILT_TAP2_MAXV - WIENER_FILT_TAP2_MINV + 1,
      WIENER_FILT_TAP2_SUBEXP_K,
      ref_wiener_info->hfilter[2] - WIENER_FILT_TAP2_MINV,
      wiener_info->hfilter[2] - WIENER_FILT_TAP2_MINV);
  memcpy(ref_wiener_info, wiener_info, sizeof(*wiener_info));
}

static void write_sgrproj_filter(const SgrprojInfo *sgrproj_info,
                                 SgrprojInfo *ref_sgrproj_info,
                                 aom_writer *wb) {
  aom_write_literal(wb, sgrproj_info->ep, SGRPROJ_PARAMS_BITS);
  const sgr_params_type *params = &av1_sgr_params[sgrproj_info->ep];

  if (params->r[0] == 0) {
    assert(sgrproj_info->xqd[0] == 0);
    aom_write_primitive_refsubexpfin(
        wb, SGRPROJ_PRJ_MAX1 - SGRPROJ_PRJ_MIN1 + 1, SGRPROJ_PRJ_SUBEXP_K,
        ref_sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1,
        sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1);
  } else if (params->r[1] == 0) {
    aom_write_primitive_refsubexpfin(
        wb, SGRPROJ_PRJ_MAX0 - SGRPROJ_PRJ_MIN0 + 1, SGRPROJ_PRJ_SUBEXP_K,
        ref_sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0,
        sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0);
  } else {
    aom_write_primitive_refsubexpfin(
        wb, SGRPROJ_PRJ_MAX0 - SGRPROJ_PRJ_MIN0 + 1, SGRPROJ_PRJ_SUBEXP_K,
        ref_sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0,
        sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0);
    aom_write_primitive_refsubexpfin(
        wb, SGRPROJ_PRJ_MAX1 - SGRPROJ_PRJ_MIN1 + 1, SGRPROJ_PRJ_SUBEXP_K,
        ref_sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1,
        sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1);
  }

  memcpy(ref_sgrproj_info, sgrproj_info, sizeof(*sgrproj_info));
}

static void loop_restoration_write_sb_coeffs(const AV1_COMMON *const cm,
                                             MACROBLOCKD *xd,
                                             const RestorationUnitInfo *rui,
                                             aom_writer *const w, int plane,
                                             FRAME_COUNTS *counts) {
  const RestorationInfo *rsi = cm->rst_info + plane;
  RestorationType frame_rtype = rsi->frame_restoration_type;
  if (frame_rtype == RESTORE_NONE) return;

  (void)counts;
  assert(!cm->all_lossless);

  const int wiener_win = (plane > 0) ? WIENER_WIN_CHROMA : WIENER_WIN;
  WienerInfo *ref_wiener_info = &xd->wiener_info[plane];
  SgrprojInfo *ref_sgrproj_info = &xd->sgrproj_info[plane];
  RestorationType unit_rtype = rui->restoration_type;

  if (frame_rtype == RESTORE_SWITCHABLE) {
    aom_write_symbol(w, unit_rtype, xd->tile_ctx->switchable_restore_cdf,
                     RESTORE_SWITCHABLE_TYPES);
#if CONFIG_ENTROPY_STATS
    ++counts->switchable_restore[unit_rtype];
#endif
    switch (unit_rtype) {
      case RESTORE_WIENER:
        write_wiener_filter(wiener_win, &rui->wiener_info, ref_wiener_info, w);
        break;
      case RESTORE_SGRPROJ:
        write_sgrproj_filter(&rui->sgrproj_info, ref_sgrproj_info, w);
        break;
      default: assert(unit_rtype == RESTORE_NONE); break;
    }
  } else if (frame_rtype == RESTORE_WIENER) {
    aom_write_symbol(w, unit_rtype != RESTORE_NONE,
                     xd->tile_ctx->wiener_restore_cdf, 2);
#if CONFIG_ENTROPY_STATS
    ++counts->wiener_restore[unit_rtype != RESTORE_NONE];
#endif
    if (unit_rtype != RESTORE_NONE) {
      write_wiener_filter(wiener_win, &rui->wiener_info, ref_wiener_info, w);
    }
  } else if (frame_rtype == RESTORE_SGRPROJ) {
    aom_write_symbol(w, unit_rtype != RESTORE_NONE,
                     xd->tile_ctx->sgrproj_restore_cdf, 2);
#if CONFIG_ENTROPY_STATS
    ++counts->sgrproj_restore[unit_rtype != RESTORE_NONE];
#endif
    if (unit_rtype != RESTORE_NONE) {
      write_sgrproj_filter(&rui->sgrproj_info, ref_sgrproj_info, w);
    }
  }
}

static void encode_loopfilter(AV1_COMMON *cm, struct aom_write_bit_buffer *wb) {
  assert(!cm->coded_lossless);
  if (cm->allow_intrabc) return;
  const int num_planes = av1_num_planes(cm);
  int i;
  struct loopfilter *lf = &cm->lf;

  // Encode the loop filter level and type
  aom_wb_write_literal(wb, lf->filter_level[0], 6);
  aom_wb_write_literal(wb, lf->filter_level[1], 6);
  if (num_planes > 1) {
    if (lf->filter_level[0] || lf->filter_level[1]) {
      aom_wb_write_literal(wb, lf->filter_level_u, 6);
      aom_wb_write_literal(wb, lf->filter_level_v, 6);
    }
  }
  aom_wb_write_literal(wb, lf->sharpness_level, 3);

  // Write out loop filter deltas applied at the MB level based on mode or
  // ref frame (if they are enabled).
  aom_wb_write_bit(wb, lf->mode_ref_delta_enabled);

  if (lf->mode_ref_delta_enabled) {
    aom_wb_write_bit(wb, lf->mode_ref_delta_update);

    if (lf->mode_ref_delta_update) {
      const RefCntBuffer *buf = get_primary_ref_frame_buf(cm);
      int8_t last_ref_deltas[REF_FRAMES];
      if (buf == NULL) {
        av1_set_default_ref_deltas(last_ref_deltas);
      } else {
        memcpy(last_ref_deltas, buf->ref_deltas, REF_FRAMES);
      }
      for (i = 0; i < REF_FRAMES; i++) {
        const int delta = lf->ref_deltas[i];
        const int changed = delta != last_ref_deltas[i];
        aom_wb_write_bit(wb, changed);
        if (changed) aom_wb_write_inv_signed_literal(wb, delta, 6);
      }

      int8_t last_mode_deltas[MAX_MODE_LF_DELTAS];
      if (buf == NULL) {
        av1_set_default_mode_deltas(last_mode_deltas);
      } else {
        memcpy(last_mode_deltas, buf->mode_deltas, MAX_MODE_LF_DELTAS);
      }
      for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        const int delta = lf->mode_deltas[i];
        const int changed = delta != last_mode_deltas[i];
        aom_wb_write_bit(wb, changed);
        if (changed) aom_wb_write_inv_signed_literal(wb, delta, 6);
      }
    }
  }
}

static void encode_cdef(const AV1_COMMON *cm, struct aom_write_bit_buffer *wb) {
  assert(!cm->coded_lossless);
  if (!cm->seq_params.enable_cdef) return;
  if (cm->allow_intrabc) return;
  const int num_planes = av1_num_planes(cm);
  int i;
  aom_wb_write_literal(wb, cm->cdef_info.cdef_damping - 3, 2);
  aom_wb_write_literal(wb, cm->cdef_info.cdef_bits, 2);
  for (i = 0; i < cm->cdef_info.nb_cdef_strengths; i++) {
    aom_wb_write_literal(wb, cm->cdef_info.cdef_strengths[i],
                         CDEF_STRENGTH_BITS);
    if (num_planes > 1)
      aom_wb_write_literal(wb, cm->cdef_info.cdef_uv_strengths[i],
                           CDEF_STRENGTH_BITS);
  }
}

static void write_delta_q(struct aom_write_bit_buffer *wb, int delta_q) {
  if (delta_q != 0) {
    aom_wb_write_bit(wb, 1);
    aom_wb_write_inv_signed_literal(wb, delta_q, 6);
  } else {
    aom_wb_write_bit(wb, 0);
  }
}

static void encode_quantization(const AV1_COMMON *const cm,
                                struct aom_write_bit_buffer *wb) {
  const int num_planes = av1_num_planes(cm);

  aom_wb_write_literal(wb, cm->base_qindex, QINDEX_BITS);
  write_delta_q(wb, cm->y_dc_delta_q);
  if (num_planes > 1) {
    int diff_uv_delta = (cm->u_dc_delta_q != cm->v_dc_delta_q) ||
                        (cm->u_ac_delta_q != cm->v_ac_delta_q);
    if (cm->seq_params.separate_uv_delta_q) aom_wb_write_bit(wb, diff_uv_delta);
    write_delta_q(wb, cm->u_dc_delta_q);
    write_delta_q(wb, cm->u_ac_delta_q);
    if (diff_uv_delta) {
      write_delta_q(wb, cm->v_dc_delta_q);
      write_delta_q(wb, cm->v_ac_delta_q);
    }
  }
  aom_wb_write_bit(wb, cm->using_qmatrix);
  if (cm->using_qmatrix) {
    aom_wb_write_literal(wb, cm->qm_y, QM_LEVEL_BITS);
    aom_wb_write_literal(wb, cm->qm_u, QM_LEVEL_BITS);
    if (!cm->seq_params.separate_uv_delta_q)
      assert(cm->qm_u == cm->qm_v);
    else
      aom_wb_write_literal(wb, cm->qm_v, QM_LEVEL_BITS);
  }
}

static void encode_segmentation(AV1_COMMON *cm, MACROBLOCKD *xd,
                                struct aom_write_bit_buffer *wb) {
  int i, j;
  struct segmentation *seg = &cm->seg;

  aom_wb_write_bit(wb, seg->enabled);
  if (!seg->enabled) return;

  // Write update flags
  if (cm->primary_ref_frame == PRIMARY_REF_NONE) {
    assert(seg->update_map == 1);
    seg->temporal_update = 0;
    assert(seg->update_data == 1);
  } else {
    aom_wb_write_bit(wb, seg->update_map);
    if (seg->update_map) {
      // Select the coding strategy (temporal or spatial)
      av1_choose_segmap_coding_method(cm, xd);
      aom_wb_write_bit(wb, seg->temporal_update);
    }
    aom_wb_write_bit(wb, seg->update_data);
  }

  // Segmentation data
  if (seg->update_data) {
    for (i = 0; i < MAX_SEGMENTS; i++) {
      for (j = 0; j < SEG_LVL_MAX; j++) {
        const int active = segfeature_active(seg, i, j);
        aom_wb_write_bit(wb, active);
        if (active) {
          const int data_max = av1_seg_feature_data_max(j);
          const int data_min = -data_max;
          const int ubits = get_unsigned_bits(data_max);
          const int data = clamp(get_segdata(seg, i, j), data_min, data_max);

          if (av1_is_segfeature_signed(j)) {
            aom_wb_write_inv_signed_literal(wb, data, ubits);
          } else {
            aom_wb_write_literal(wb, data, ubits);
          }
        }
      }
    }
  }
}

static void write_frame_interp_filter(InterpFilter filter,
                                      struct aom_write_bit_buffer *wb) {
  aom_wb_write_bit(wb, filter == SWITCHABLE);
  if (filter != SWITCHABLE)
    aom_wb_write_literal(wb, filter, LOG_SWITCHABLE_FILTERS);
}

// Same function as write_uniform but writing to uncompresses header wb
static void wb_write_uniform(struct aom_write_bit_buffer *wb, int n, int v) {
  const int l = get_unsigned_bits(n);
  const int m = (1 << l) - n;
  if (l == 0) return;
  if (v < m) {
    aom_wb_write_literal(wb, v, l - 1);
  } else {
    aom_wb_write_literal(wb, m + ((v - m) >> 1), l - 1);
    aom_wb_write_literal(wb, (v - m) & 1, 1);
  }
}

static void write_tile_info_max_tile(const AV1_COMMON *const cm,
                                     struct aom_write_bit_buffer *wb) {
  int width_mi = ALIGN_POWER_OF_TWO(cm->mi_cols, cm->seq_params.mib_size_log2);
  int height_mi = ALIGN_POWER_OF_TWO(cm->mi_rows, cm->seq_params.mib_size_log2);
  int width_sb = width_mi >> cm->seq_params.mib_size_log2;
  int height_sb = height_mi >> cm->seq_params.mib_size_log2;
  int size_sb, i;

  aom_wb_write_bit(wb, cm->uniform_tile_spacing_flag);

  if (cm->uniform_tile_spacing_flag) {
    // Uniform spaced tiles with power-of-two number of rows and columns
    // tile columns
    int ones = cm->log2_tile_cols - cm->min_log2_tile_cols;
    while (ones--) {
      aom_wb_write_bit(wb, 1);
    }
    if (cm->log2_tile_cols < cm->max_log2_tile_cols) {
      aom_wb_write_bit(wb, 0);
    }

    // rows
    ones = cm->log2_tile_rows - cm->min_log2_tile_rows;
    while (ones--) {
      aom_wb_write_bit(wb, 1);
    }
    if (cm->log2_tile_rows < cm->max_log2_tile_rows) {
      aom_wb_write_bit(wb, 0);
    }
  } else {
    // Explicit tiles with configurable tile widths and heights
    // columns
    for (i = 0; i < cm->tile_cols; i++) {
      size_sb = cm->tile_col_start_sb[i + 1] - cm->tile_col_start_sb[i];
      wb_write_uniform(wb, AOMMIN(width_sb, cm->max_tile_width_sb),
                       size_sb - 1);
      width_sb -= size_sb;
    }
    assert(width_sb == 0);

    // rows
    for (i = 0; i < cm->tile_rows; i++) {
      size_sb = cm->tile_row_start_sb[i + 1] - cm->tile_row_start_sb[i];
      wb_write_uniform(wb, AOMMIN(height_sb, cm->max_tile_height_sb),
                       size_sb - 1);
      height_sb -= size_sb;
    }
    assert(height_sb == 0);
  }
}

static void write_tile_info(const AV1_COMMON *const cm,
                            struct aom_write_bit_buffer *saved_wb,
                            struct aom_write_bit_buffer *wb) {
  write_tile_info_max_tile(cm, wb);

  *saved_wb = *wb;
  if (cm->tile_rows * cm->tile_cols > 1) {
    // tile id used for cdf update
    aom_wb_write_literal(wb, 0, cm->log2_tile_cols + cm->log2_tile_rows);
    // Number of bytes in tile size - 1
    aom_wb_write_literal(wb, 3, 2);
  }
}

static void write_ext_tile_info(const AV1_COMMON *const cm,
                                struct aom_write_bit_buffer *saved_wb,
                                struct aom_write_bit_buffer *wb) {
  // This information is stored as a separate byte.
  int mod = wb->bit_offset % CHAR_BIT;
  if (mod > 0) aom_wb_write_literal(wb, 0, CHAR_BIT - mod);
  assert(aom_wb_is_byte_aligned(wb));

  *saved_wb = *wb;
  if (cm->tile_rows * cm->tile_cols > 1) {
    // Note that the last item in the uncompressed header is the data
    // describing tile configuration.
    // Number of bytes in tile column size - 1
    aom_wb_write_literal(wb, 0, 2);
    // Number of bytes in tile size - 1
    aom_wb_write_literal(wb, 0, 2);
  }
}

// Stores the location and size of a tile's data in the bitstream.  Used for
// later identifying identical tiles
typedef struct TileBufferEnc {
  uint8_t *data;
  size_t size;
} TileBufferEnc;

static INLINE int find_identical_tile(
    const int tile_row, const int tile_col,
    TileBufferEnc (*const tile_buffers)[MAX_TILE_COLS]) {
  const MV32 candidate_offset[1] = { { 1, 0 } };
  const uint8_t *const cur_tile_data =
      tile_buffers[tile_row][tile_col].data + 4;
  const size_t cur_tile_size = tile_buffers[tile_row][tile_col].size;

  int i;

  if (tile_row == 0) return 0;

  // (TODO: yunqingwang) For now, only above tile is checked and used.
  // More candidates such as left tile can be added later.
  for (i = 0; i < 1; i++) {
    int row_offset = candidate_offset[0].row;
    int col_offset = candidate_offset[0].col;
    int row = tile_row - row_offset;
    int col = tile_col - col_offset;
    const uint8_t *tile_data;
    TileBufferEnc *candidate;

    if (row < 0 || col < 0) continue;

    const uint32_t tile_hdr = mem_get_le32(tile_buffers[row][col].data);

    // Read out tile-copy-mode bit:
    if ((tile_hdr >> 31) == 1) {
      // The candidate is a copy tile itself: the offset is stored in bits
      // 30 through 24 inclusive.
      row_offset += (tile_hdr >> 24) & 0x7f;
      row = tile_row - row_offset;
    }

    candidate = &tile_buffers[row][col];

    if (row_offset >= 128 || candidate->size != cur_tile_size) continue;

    tile_data = candidate->data + 4;

    if (memcmp(tile_data, cur_tile_data, cur_tile_size) != 0) continue;

    // Identical tile found
    assert(row_offset > 0);
    return row_offset;
  }

  // No identical tile found
  return 0;
}

static void write_render_size(const AV1_COMMON *cm,
                              struct aom_write_bit_buffer *wb) {
  const int scaling_active = av1_resize_scaled(cm);
  aom_wb_write_bit(wb, scaling_active);
  if (scaling_active) {
    aom_wb_write_literal(wb, cm->render_width - 1, 16);
    aom_wb_write_literal(wb, cm->render_height - 1, 16);
  }
}

static void write_superres_scale(const AV1_COMMON *const cm,
                                 struct aom_write_bit_buffer *wb) {
  const SequenceHeader *const seq_params = &cm->seq_params;
  if (!seq_params->enable_superres) {
    assert(cm->superres_scale_denominator == SCALE_NUMERATOR);
    return;
  }

  // First bit is whether to to scale or not
  if (cm->superres_scale_denominator == SCALE_NUMERATOR) {
    aom_wb_write_bit(wb, 0);  // no scaling
  } else {
    aom_wb_write_bit(wb, 1);  // scaling, write scale factor
    assert(cm->superres_scale_denominator >= SUPERRES_SCALE_DENOMINATOR_MIN);
    assert(cm->superres_scale_denominator <
           SUPERRES_SCALE_DENOMINATOR_MIN + (1 << SUPERRES_SCALE_BITS));
    aom_wb_write_literal(
        wb, cm->superres_scale_denominator - SUPERRES_SCALE_DENOMINATOR_MIN,
        SUPERRES_SCALE_BITS);
  }
}

static void write_frame_size(const AV1_COMMON *cm, int frame_size_override,
                             struct aom_write_bit_buffer *wb) {
  const int coded_width = cm->superres_upscaled_width - 1;
  const int coded_height = cm->superres_upscaled_height - 1;

  if (frame_size_override) {
    const SequenceHeader *seq_params = &cm->seq_params;
    int num_bits_width = seq_params->num_bits_width;
    int num_bits_height = seq_params->num_bits_height;
    aom_wb_write_literal(wb, coded_width, num_bits_width);
    aom_wb_write_literal(wb, coded_height, num_bits_height);
  }

  write_superres_scale(cm, wb);
  write_render_size(cm, wb);
}

static void write_frame_size_with_refs(const AV1_COMMON *const cm,
                                       struct aom_write_bit_buffer *wb) {
  int found = 0;

  MV_REFERENCE_FRAME ref_frame;
  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    const YV12_BUFFER_CONFIG *cfg = get_ref_frame_yv12_buf(cm, ref_frame);

    if (cfg != NULL) {
      found = cm->superres_upscaled_width == cfg->y_crop_width &&
              cm->superres_upscaled_height == cfg->y_crop_height;
      found &= cm->render_width == cfg->render_width &&
               cm->render_height == cfg->render_height;
    }
    aom_wb_write_bit(wb, found);
    if (found) {
      write_superres_scale(cm, wb);
      break;
    }
  }

  if (!found) {
    int frame_size_override = 1;  // Always equal to 1 in this function
    write_frame_size(cm, frame_size_override, wb);
  }
}

static void write_profile(BITSTREAM_PROFILE profile,
                          struct aom_write_bit_buffer *wb) {
  assert(profile >= PROFILE_0 && profile < MAX_PROFILES);
  aom_wb_write_literal(wb, profile, PROFILE_BITS);
}

static void write_bitdepth(const SequenceHeader *const seq_params,
                           struct aom_write_bit_buffer *wb) {
  // Profile 0/1: [0] for 8 bit, [1]  10-bit
  // Profile   2: [0] for 8 bit, [10] 10-bit, [11] - 12-bit
  aom_wb_write_bit(wb, seq_params->bit_depth == AOM_BITS_8 ? 0 : 1);
  if (seq_params->profile == PROFILE_2 && seq_params->bit_depth != AOM_BITS_8) {
    aom_wb_write_bit(wb, seq_params->bit_depth == AOM_BITS_10 ? 0 : 1);
  }
}

static void write_color_config(const SequenceHeader *const seq_params,
                               struct aom_write_bit_buffer *wb) {
  write_bitdepth(seq_params, wb);
  const int is_monochrome = seq_params->monochrome;
  // monochrome bit
  if (seq_params->profile != PROFILE_1)
    aom_wb_write_bit(wb, is_monochrome);
  else
    assert(!is_monochrome);
  if (seq_params->color_primaries == AOM_CICP_CP_UNSPECIFIED &&
      seq_params->transfer_characteristics == AOM_CICP_TC_UNSPECIFIED &&
      seq_params->matrix_coefficients == AOM_CICP_MC_UNSPECIFIED) {
    aom_wb_write_bit(wb, 0);  // No color description present
  } else {
    aom_wb_write_bit(wb, 1);  // Color description present
    aom_wb_write_literal(wb, seq_params->color_primaries, 8);
    aom_wb_write_literal(wb, seq_params->transfer_characteristics, 8);
    aom_wb_write_literal(wb, seq_params->matrix_coefficients, 8);
  }
  if (is_monochrome) {
    // 0: [16, 235] (i.e. xvYCC), 1: [0, 255]
    aom_wb_write_bit(wb, seq_params->color_range);
    return;
  }
  if (seq_params->color_primaries == AOM_CICP_CP_BT_709 &&
      seq_params->transfer_characteristics == AOM_CICP_TC_SRGB &&
      seq_params->matrix_coefficients == AOM_CICP_MC_IDENTITY) {
    assert(seq_params->subsampling_x == 0 && seq_params->subsampling_y == 0);
    assert(seq_params->profile == PROFILE_1 ||
           (seq_params->profile == PROFILE_2 &&
            seq_params->bit_depth == AOM_BITS_12));
  } else {
    // 0: [16, 235] (i.e. xvYCC), 1: [0, 255]
    aom_wb_write_bit(wb, seq_params->color_range);
    if (seq_params->profile == PROFILE_0) {
      // 420 only
      assert(seq_params->subsampling_x == 1 && seq_params->subsampling_y == 1);
    } else if (seq_params->profile == PROFILE_1) {
      // 444 only
      assert(seq_params->subsampling_x == 0 && seq_params->subsampling_y == 0);
    } else if (seq_params->profile == PROFILE_2) {
      if (seq_params->bit_depth == AOM_BITS_12) {
        // 420, 444 or 422
        aom_wb_write_bit(wb, seq_params->subsampling_x);
        if (seq_params->subsampling_x == 0) {
          assert(seq_params->subsampling_y == 0 &&
                 "4:4:0 subsampling not allowed in AV1");
        } else {
          aom_wb_write_bit(wb, seq_params->subsampling_y);
        }
      } else {
        // 422 only
        assert(seq_params->subsampling_x == 1 &&
               seq_params->subsampling_y == 0);
      }
    }
    if (seq_params->matrix_coefficients == AOM_CICP_MC_IDENTITY) {
      assert(seq_params->subsampling_x == 0 && seq_params->subsampling_y == 0);
    }
    if (seq_params->subsampling_x == 1 && seq_params->subsampling_y == 1) {
      aom_wb_write_literal(wb, seq_params->chroma_sample_position, 2);
    }
  }
  aom_wb_write_bit(wb, seq_params->separate_uv_delta_q);
}

static void write_timing_info_header(AV1_COMMON *const cm,
                                     struct aom_write_bit_buffer *wb) {
  aom_wb_write_unsigned_literal(wb, cm->timing_info.num_units_in_display_tick,
                                32);  // Number of units in tick
  aom_wb_write_unsigned_literal(wb, cm->timing_info.time_scale,
                                32);  // Time scale
  aom_wb_write_bit(
      wb,
      cm->timing_info.equal_picture_interval);  // Equal picture interval bit
  if (cm->timing_info.equal_picture_interval) {
    aom_wb_write_uvlc(
        wb,
        cm->timing_info.num_ticks_per_picture - 1);  // ticks per picture
  }
}

static void write_decoder_model_info(AV1_COMMON *const cm,
                                     struct aom_write_bit_buffer *wb) {
  aom_wb_write_literal(
      wb, cm->buffer_model.encoder_decoder_buffer_delay_length - 1, 5);
  aom_wb_write_unsigned_literal(wb, cm->buffer_model.num_units_in_decoding_tick,
                                32);  // Number of units in decoding tick
  aom_wb_write_literal(wb, cm->buffer_model.buffer_removal_time_length - 1, 5);
  aom_wb_write_literal(wb, cm->buffer_model.frame_presentation_time_length - 1,
                       5);
}

static void write_dec_model_op_parameters(AV1_COMMON *const cm,
                                          struct aom_write_bit_buffer *wb,
                                          int op_num) {
  if (op_num > MAX_NUM_OPERATING_POINTS)
    aom_internal_error(
        &cm->error, AOM_CODEC_UNSUP_BITSTREAM,
        "Encoder does not support %d decoder model operating points", op_num);

  //  aom_wb_write_bit(wb, cm->op_params[op_num].has_parameters);
  //  if (!cm->op_params[op_num].has_parameters) return;

  aom_wb_write_unsigned_literal(
      wb, cm->op_params[op_num].decoder_buffer_delay,
      cm->buffer_model.encoder_decoder_buffer_delay_length);

  aom_wb_write_unsigned_literal(
      wb, cm->op_params[op_num].encoder_buffer_delay,
      cm->buffer_model.encoder_decoder_buffer_delay_length);

  aom_wb_write_bit(wb, cm->op_params[op_num].low_delay_mode_flag);

  cm->op_frame_timing[op_num].buffer_removal_time =
      0;  // reset the decoded frame counter
}

static void write_tu_pts_info(AV1_COMMON *const cm,
                              struct aom_write_bit_buffer *wb) {
  aom_wb_write_unsigned_literal(
      wb, cm->frame_presentation_time,
      cm->buffer_model.frame_presentation_time_length);
}

static void write_film_grain_params(const AV1_COMP *const cpi,
                                    struct aom_write_bit_buffer *wb) {
  const AV1_COMMON *const cm = &cpi->common;
  const aom_film_grain_t *const pars = &cm->cur_frame->film_grain_params;

  aom_wb_write_bit(wb, pars->apply_grain);
  if (!pars->apply_grain) return;

  aom_wb_write_literal(wb, pars->random_seed, 16);

  if (cm->current_frame.frame_type == INTER_FRAME)
    aom_wb_write_bit(wb, pars->update_parameters);

  if (!pars->update_parameters) {
    int ref_frame, ref_idx;
    for (ref_frame = LAST_FRAME; ref_frame < REF_FRAMES; ref_frame++) {
      ref_idx = get_ref_frame_map_idx(cm, ref_frame);
      assert(ref_idx != INVALID_IDX);
      const RefCntBuffer *const buf = cm->ref_frame_map[ref_idx];
      if (buf->film_grain_params_present &&
          av1_check_grain_params_equiv(pars, &buf->film_grain_params)) {
        break;
      }
    }
    assert(ref_frame < REF_FRAMES);
    aom_wb_write_literal(wb, ref_idx, 3);
    return;
  }

  // Scaling functions parameters
  aom_wb_write_literal(wb, pars->num_y_points, 4);  // max 14
  for (int i = 0; i < pars->num_y_points; i++) {
    aom_wb_write_literal(wb, pars->scaling_points_y[i][0], 8);
    aom_wb_write_literal(wb, pars->scaling_points_y[i][1], 8);
  }

  if (!cm->seq_params.monochrome) {
    aom_wb_write_bit(wb, pars->chroma_scaling_from_luma);
  } else {
    assert(!pars->chroma_scaling_from_luma);
  }

  if (cm->seq_params.monochrome || pars->chroma_scaling_from_luma ||
      ((cm->seq_params.subsampling_x == 1) &&
       (cm->seq_params.subsampling_y == 1) && (pars->num_y_points == 0))) {
    assert(pars->num_cb_points == 0 && pars->num_cr_points == 0);
  } else {
    aom_wb_write_literal(wb, pars->num_cb_points, 4);  // max 10
    for (int i = 0; i < pars->num_cb_points; i++) {
      aom_wb_write_literal(wb, pars->scaling_points_cb[i][0], 8);
      aom_wb_write_literal(wb, pars->scaling_points_cb[i][1], 8);
    }

    aom_wb_write_literal(wb, pars->num_cr_points, 4);  // max 10
    for (int i = 0; i < pars->num_cr_points; i++) {
      aom_wb_write_literal(wb, pars->scaling_points_cr[i][0], 8);
      aom_wb_write_literal(wb, pars->scaling_points_cr[i][1], 8);
    }
  }

  aom_wb_write_literal(wb, pars->scaling_shift - 8, 2);  // 8 + value

  // AR coefficients
  // Only sent if the corresponsing scaling function has
  // more than 0 points

  aom_wb_write_literal(wb, pars->ar_coeff_lag, 2);

  int num_pos_luma = 2 * pars->ar_coeff_lag * (pars->ar_coeff_lag + 1);
  int num_pos_chroma = num_pos_luma;
  if (pars->num_y_points > 0) ++num_pos_chroma;

  if (pars->num_y_points)
    for (int i = 0; i < num_pos_luma; i++)
      aom_wb_write_literal(wb, pars->ar_coeffs_y[i] + 128, 8);

  if (pars->num_cb_points || pars->chroma_scaling_from_luma)
    for (int i = 0; i < num_pos_chroma; i++)
      aom_wb_write_literal(wb, pars->ar_coeffs_cb[i] + 128, 8);

  if (pars->num_cr_points || pars->chroma_scaling_from_luma)
    for (int i = 0; i < num_pos_chroma; i++)
      aom_wb_write_literal(wb, pars->ar_coeffs_cr[i] + 128, 8);

  aom_wb_write_literal(wb, pars->ar_coeff_shift - 6, 2);  // 8 + value

  aom_wb_write_literal(wb, pars->grain_scale_shift, 2);

  if (pars->num_cb_points) {
    aom_wb_write_literal(wb, pars->cb_mult, 8);
    aom_wb_write_literal(wb, pars->cb_luma_mult, 8);
    aom_wb_write_literal(wb, pars->cb_offset, 9);
  }

  if (pars->num_cr_points) {
    aom_wb_write_literal(wb, pars->cr_mult, 8);
    aom_wb_write_literal(wb, pars->cr_luma_mult, 8);
    aom_wb_write_literal(wb, pars->cr_offset, 9);
  }

  aom_wb_write_bit(wb, pars->overlap_flag);

  aom_wb_write_bit(wb, pars->clip_to_restricted_range);
}

static void write_sb_size(const SequenceHeader *const seq_params,
                          struct aom_write_bit_buffer *wb) {
  (void)seq_params;
  (void)wb;
  assert(seq_params->mib_size == mi_size_wide[seq_params->sb_size]);
  assert(seq_params->mib_size == 1 << seq_params->mib_size_log2);
  assert(seq_params->sb_size == BLOCK_128X128 ||
         seq_params->sb_size == BLOCK_64X64);
  aom_wb_write_bit(wb, seq_params->sb_size == BLOCK_128X128 ? 1 : 0);
}

static void write_sequence_header(const SequenceHeader *const seq_params,
                                  struct aom_write_bit_buffer *wb) {
  aom_wb_write_literal(wb, seq_params->num_bits_width - 1, 4);
  aom_wb_write_literal(wb, seq_params->num_bits_height - 1, 4);
  aom_wb_write_literal(wb, seq_params->max_frame_width - 1,
                       seq_params->num_bits_width);
  aom_wb_write_literal(wb, seq_params->max_frame_height - 1,
                       seq_params->num_bits_height);

  if (!seq_params->reduced_still_picture_hdr) {
    aom_wb_write_bit(wb, seq_params->frame_id_numbers_present_flag);
    if (seq_params->frame_id_numbers_present_flag) {
      // We must always have delta_frame_id_length < frame_id_length,
      // in order for a frame to be referenced with a unique delta.
      // Avoid wasting bits by using a coding that enforces this restriction.
      aom_wb_write_literal(wb, seq_params->delta_frame_id_length - 2, 4);
      aom_wb_write_literal(
          wb,
          seq_params->frame_id_length - seq_params->delta_frame_id_length - 1,
          3);
    }
  }

  write_sb_size(seq_params, wb);

  aom_wb_write_bit(wb, seq_params->enable_filter_intra);
  aom_wb_write_bit(wb, seq_params->enable_intra_edge_filter);

  if (!seq_params->reduced_still_picture_hdr) {
    aom_wb_write_bit(wb, seq_params->enable_interintra_compound);
    aom_wb_write_bit(wb, seq_params->enable_masked_compound);
    aom_wb_write_bit(wb, seq_params->enable_warped_motion);
    aom_wb_write_bit(wb, seq_params->enable_dual_filter);

    aom_wb_write_bit(wb, seq_params->order_hint_info.enable_order_hint);

    if (seq_params->order_hint_info.enable_order_hint) {
      aom_wb_write_bit(wb, seq_params->order_hint_info.enable_dist_wtd_comp);
      aom_wb_write_bit(wb, seq_params->order_hint_info.enable_ref_frame_mvs);
    }
    if (seq_params->force_screen_content_tools == 2) {
      aom_wb_write_bit(wb, 1);
    } else {
      aom_wb_write_bit(wb, 0);
      aom_wb_write_bit(wb, seq_params->force_screen_content_tools);
    }
    if (seq_params->force_screen_content_tools > 0) {
      if (seq_params->force_integer_mv == 2) {
        aom_wb_write_bit(wb, 1);
      } else {
        aom_wb_write_bit(wb, 0);
        aom_wb_write_bit(wb, seq_params->force_integer_mv);
      }
    } else {
      assert(seq_params->force_integer_mv == 2);
    }
    if (seq_params->order_hint_info.enable_order_hint)
      aom_wb_write_literal(
          wb, seq_params->order_hint_info.order_hint_bits_minus_1, 3);
  }

  aom_wb_write_bit(wb, seq_params->enable_superres);
  aom_wb_write_bit(wb, seq_params->enable_cdef);
  aom_wb_write_bit(wb, seq_params->enable_restoration);
}

static void write_global_motion_params(const WarpedMotionParams *params,
                                       const WarpedMotionParams *ref_params,
                                       struct aom_write_bit_buffer *wb,
                                       int allow_hp) {
  const TransformationType type = params->wmtype;

  aom_wb_write_bit(wb, type != IDENTITY);
  if (type != IDENTITY) {
    aom_wb_write_bit(wb, type == ROTZOOM);
    if (type != ROTZOOM) aom_wb_write_bit(wb, type == TRANSLATION);
  }

  if (type >= ROTZOOM) {
    aom_wb_write_signed_primitive_refsubexpfin(
        wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
        (ref_params->wmmat[2] >> GM_ALPHA_PREC_DIFF) -
            (1 << GM_ALPHA_PREC_BITS),
        (params->wmmat[2] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));
    aom_wb_write_signed_primitive_refsubexpfin(
        wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
        (ref_params->wmmat[3] >> GM_ALPHA_PREC_DIFF),
        (params->wmmat[3] >> GM_ALPHA_PREC_DIFF));
  }

  if (type >= AFFINE) {
    aom_wb_write_signed_primitive_refsubexpfin(
        wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
        (ref_params->wmmat[4] >> GM_ALPHA_PREC_DIFF),
        (params->wmmat[4] >> GM_ALPHA_PREC_DIFF));
    aom_wb_write_signed_primitive_refsubexpfin(
        wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
        (ref_params->wmmat[5] >> GM_ALPHA_PREC_DIFF) -
            (1 << GM_ALPHA_PREC_BITS),
        (params->wmmat[5] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));
  }

  if (type >= TRANSLATION) {
    const int trans_bits = (type == TRANSLATION)
                               ? GM_ABS_TRANS_ONLY_BITS - !allow_hp
                               : GM_ABS_TRANS_BITS;
    const int trans_prec_diff = (type == TRANSLATION)
                                    ? GM_TRANS_ONLY_PREC_DIFF + !allow_hp
                                    : GM_TRANS_PREC_DIFF;
    aom_wb_write_signed_primitive_refsubexpfin(
        wb, (1 << trans_bits) + 1, SUBEXPFIN_K,
        (ref_params->wmmat[0] >> trans_prec_diff),
        (params->wmmat[0] >> trans_prec_diff));
    aom_wb_write_signed_primitive_refsubexpfin(
        wb, (1 << trans_bits) + 1, SUBEXPFIN_K,
        (ref_params->wmmat[1] >> trans_prec_diff),
        (params->wmmat[1] >> trans_prec_diff));
  }
}

static void write_global_motion(AV1_COMP *cpi,
                                struct aom_write_bit_buffer *wb) {
  AV1_COMMON *const cm = &cpi->common;
  int frame;
  for (frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame) {
    const WarpedMotionParams *ref_params =
        cm->prev_frame ? &cm->prev_frame->global_motion[frame]
                       : &default_warp_params;
    write_global_motion_params(&cm->global_motion[frame], ref_params, wb,
                               cm->allow_high_precision_mv);
    // TODO(sarahparker, debargha): The logic in the commented out code below
    // does not work currently and causes mismatches when resize is on.
    // Fix it before turning the optimization back on.
    /*
    YV12_BUFFER_CONFIG *ref_buf = get_ref_frame_yv12_buf(cpi, frame);
    if (cpi->source->y_crop_width == ref_buf->y_crop_width &&
        cpi->source->y_crop_height == ref_buf->y_crop_height) {
      write_global_motion_params(&cm->global_motion[frame],
                                 &cm->prev_frame->global_motion[frame], wb,
                                 cm->allow_high_precision_mv);
    } else {
      assert(cm->global_motion[frame].wmtype == IDENTITY &&
             "Invalid warp type for frames of different resolutions");
    }
    */
    /*
    printf("Frame %d/%d: Enc Ref %d: %d %d %d %d\n",
           cm->current_frame.frame_number, cm->show_frame, frame,
           cm->global_motion[frame].wmmat[0],
           cm->global_motion[frame].wmmat[1], cm->global_motion[frame].wmmat[2],
           cm->global_motion[frame].wmmat[3]);
           */
  }
}

static int check_frame_refs_short_signaling(AV1_COMMON *const cm) {
  // Check whether all references are distinct frames.
  const RefCntBuffer *seen_bufs[FRAME_BUFFERS] = { NULL };
  int num_refs = 0;
  for (int ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    const RefCntBuffer *const buf = get_ref_frame_buf(cm, ref_frame);
    if (buf != NULL) {
      int seen = 0;
      for (int i = 0; i < num_refs; i++) {
        if (seen_bufs[i] == buf) {
          seen = 1;
          break;
        }
      }
      if (!seen) seen_bufs[num_refs++] = buf;
    }
  }

  // We only turn on frame_refs_short_signaling when all references are
  // distinct.
  if (num_refs < INTER_REFS_PER_FRAME) {
    // It indicates that there exist more than one reference frame pointing to
    // the same reference buffer, i.e. two or more references are duplicate.
    return 0;
  }

  // Check whether the encoder side ref frame choices are aligned with that to
  // be derived at the decoder side.
  int remapped_ref_idx_decoder[REF_FRAMES];

  const int lst_map_idx = get_ref_frame_map_idx(cm, LAST_FRAME);
  const int gld_map_idx = get_ref_frame_map_idx(cm, GOLDEN_FRAME);

  // Set up the frame refs mapping indexes according to the
  // frame_refs_short_signaling policy.
  av1_set_frame_refs(cm, remapped_ref_idx_decoder, lst_map_idx, gld_map_idx);

  // We only turn on frame_refs_short_signaling when the encoder side decision
  // on ref frames is identical to that at the decoder side.
  int frame_refs_short_signaling = 1;
  for (int ref_idx = 0; ref_idx < INTER_REFS_PER_FRAME; ++ref_idx) {
    // Compare the buffer index between two reference frames indexed
    // respectively by the encoder and the decoder side decisions.
    RefCntBuffer *ref_frame_buf_new = NULL;
    if (remapped_ref_idx_decoder[ref_idx] != INVALID_IDX) {
      ref_frame_buf_new = cm->ref_frame_map[remapped_ref_idx_decoder[ref_idx]];
    }
    if (get_ref_frame_buf(cm, LAST_FRAME + ref_idx) != ref_frame_buf_new) {
      frame_refs_short_signaling = 0;
      break;
    }
  }

#if 0   // For debug
  printf("\nFrame=%d: \n", cm->current_frame.frame_number);
  printf("***frame_refs_short_signaling=%d\n", frame_refs_short_signaling);
  for (int ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    printf("enc_ref(map_idx=%d)=%d, vs. "
        "dec_ref(map_idx=%d)=%d\n",
        get_ref_frame_map_idx(cm, ref_frame), ref_frame,
        cm->remapped_ref_idx[ref_frame - LAST_FRAME],
        ref_frame);
  }
#endif  // 0

  return frame_refs_short_signaling;
}

// New function based on HLS R18
static void write_uncompressed_header_obu(AV1_COMP *cpi,
                                          struct aom_write_bit_buffer *saved_wb,
                                          struct aom_write_bit_buffer *wb) {
  AV1_COMMON *const cm = &cpi->common;
  const SequenceHeader *const seq_params = &cm->seq_params;
  MACROBLOCKD *const xd = &cpi->td.mb.e_mbd;
  CurrentFrame *const current_frame = &cm->current_frame;

  current_frame->frame_refs_short_signaling = 0;

  if (seq_params->still_picture) {
    assert(cm->show_existing_frame == 0);
    assert(cm->show_frame == 1);
    assert(current_frame->frame_type == KEY_FRAME);
  }
  if (!seq_params->reduced_still_picture_hdr) {
    if (encode_show_existing_frame(cm)) {
      aom_wb_write_bit(wb, 1);  // show_existing_frame
      aom_wb_write_literal(wb, cpi->existing_fb_idx_to_show, 3);

      if (seq_params->decoder_model_info_present_flag &&
          cm->timing_info.equal_picture_interval == 0) {
        write_tu_pts_info(cm, wb);
      }
      if (seq_params->frame_id_numbers_present_flag) {
        int frame_id_len = seq_params->frame_id_length;
        int display_frame_id = cm->ref_frame_id[cpi->existing_fb_idx_to_show];
        aom_wb_write_literal(wb, display_frame_id, frame_id_len);
      }
      return;
    } else {
      aom_wb_write_bit(wb, 0);  // show_existing_frame
    }

    aom_wb_write_literal(wb, current_frame->frame_type, 2);

    aom_wb_write_bit(wb, cm->show_frame);
    if (cm->show_frame) {
      if (seq_params->decoder_model_info_present_flag &&
          cm->timing_info.equal_picture_interval == 0)
        write_tu_pts_info(cm, wb);
    } else {
      aom_wb_write_bit(wb, cm->showable_frame);
    }
    if (frame_is_sframe(cm)) {
      assert(cm->error_resilient_mode);
    } else if (!(current_frame->frame_type == KEY_FRAME && cm->show_frame)) {
      aom_wb_write_bit(wb, cm->error_resilient_mode);
    }
  }
  aom_wb_write_bit(wb, cm->disable_cdf_update);

  if (seq_params->force_screen_content_tools == 2) {
    aom_wb_write_bit(wb, cm->allow_screen_content_tools);
  } else {
    assert(cm->allow_screen_content_tools ==
           seq_params->force_screen_content_tools);
  }

  if (cm->allow_screen_content_tools) {
    if (seq_params->force_integer_mv == 2) {
      aom_wb_write_bit(wb, cm->cur_frame_force_integer_mv);
    } else {
      assert(cm->cur_frame_force_integer_mv == seq_params->force_integer_mv);
    }
  } else {
    assert(cm->cur_frame_force_integer_mv == 0);
  }

  int frame_size_override_flag = 0;

  if (seq_params->reduced_still_picture_hdr) {
    assert(cm->superres_upscaled_width == seq_params->max_frame_width &&
           cm->superres_upscaled_height == seq_params->max_frame_height);
  } else {
    if (seq_params->frame_id_numbers_present_flag) {
      int frame_id_len = seq_params->frame_id_length;
      aom_wb_write_literal(wb, cm->current_frame_id, frame_id_len);
    }

    if (cm->superres_upscaled_width > seq_params->max_frame_width ||
        cm->superres_upscaled_height > seq_params->max_frame_height) {
      aom_internal_error(&cm->error, AOM_CODEC_UNSUP_BITSTREAM,
                         "Frame dimensions are larger than the maximum values");
    }

    frame_size_override_flag =
        frame_is_sframe(cm)
            ? 1
            : (cm->superres_upscaled_width != seq_params->max_frame_width ||
               cm->superres_upscaled_height != seq_params->max_frame_height);
    if (!frame_is_sframe(cm)) aom_wb_write_bit(wb, frame_size_override_flag);

    if (seq_params->order_hint_info.enable_order_hint)
      aom_wb_write_literal(
          wb, current_frame->order_hint,
          seq_params->order_hint_info.order_hint_bits_minus_1 + 1);

    if (!cm->error_resilient_mode && !frame_is_intra_only(cm)) {
      aom_wb_write_literal(wb, cm->primary_ref_frame, PRIMARY_REF_BITS);
    }
  }

  if (seq_params->decoder_model_info_present_flag) {
    aom_wb_write_bit(wb, cm->buffer_removal_time_present);
    if (cm->buffer_removal_time_present) {
      for (int op_num = 0;
           op_num < seq_params->operating_points_cnt_minus_1 + 1; op_num++) {
        if (cm->op_params[op_num].decoder_model_param_present_flag) {
          if (((seq_params->operating_point_idc[op_num] >>
                cm->temporal_layer_id) &
                   0x1 &&
               (seq_params->operating_point_idc[op_num] >>
                (cm->spatial_layer_id + 8)) &
                   0x1) ||
              seq_params->operating_point_idc[op_num] == 0) {
            aom_wb_write_unsigned_literal(
                wb, cm->op_frame_timing[op_num].buffer_removal_time,
                cm->buffer_model.buffer_removal_time_length);
            cm->op_frame_timing[op_num].buffer_removal_time++;
            if (cm->op_frame_timing[op_num].buffer_removal_time == 0) {
              aom_internal_error(&cm->error, AOM_CODEC_UNSUP_BITSTREAM,
                                 "buffer_removal_time overflowed");
            }
          }
        }
      }
    }
  }

  // Shown keyframes and switch-frames automatically refreshes all reference
  // frames.  For all other frame types, we need to write refresh_frame_flags.
  if ((current_frame->frame_type == KEY_FRAME && !cm->show_frame) ||
      current_frame->frame_type == INTER_FRAME ||
      current_frame->frame_type == INTRA_ONLY_FRAME)
    aom_wb_write_literal(wb, current_frame->refresh_frame_flags, REF_FRAMES);

  if (!frame_is_intra_only(cm) || current_frame->refresh_frame_flags != 0xff) {
    // Write all ref frame order hints if error_resilient_mode == 1
    if (cm->error_resilient_mode &&
        seq_params->order_hint_info.enable_order_hint) {
      for (int ref_idx = 0; ref_idx < REF_FRAMES; ref_idx++) {
        aom_wb_write_literal(
            wb, cm->ref_frame_map[ref_idx]->order_hint,
            seq_params->order_hint_info.order_hint_bits_minus_1 + 1);
      }
    }
  }

  if (current_frame->frame_type == KEY_FRAME) {
    write_frame_size(cm, frame_size_override_flag, wb);
    assert(!av1_superres_scaled(cm) || !cm->allow_intrabc);
    if (cm->allow_screen_content_tools && !av1_superres_scaled(cm))
      aom_wb_write_bit(wb, cm->allow_intrabc);
  } else {
    if (current_frame->frame_type == INTRA_ONLY_FRAME) {
      write_frame_size(cm, frame_size_override_flag, wb);
      assert(!av1_superres_scaled(cm) || !cm->allow_intrabc);
      if (cm->allow_screen_content_tools && !av1_superres_scaled(cm))
        aom_wb_write_bit(wb, cm->allow_intrabc);
    } else if (current_frame->frame_type == INTER_FRAME ||
               frame_is_sframe(cm)) {
      MV_REFERENCE_FRAME ref_frame;

      // NOTE: Error resilient mode turns off frame_refs_short_signaling
      //       automatically.
#define FRAME_REFS_SHORT_SIGNALING 0
#if FRAME_REFS_SHORT_SIGNALING
      current_frame->frame_refs_short_signaling =
          seq_params->order_hint_info.enable_order_hint;
#endif  // FRAME_REFS_SHORT_SIGNALING

      if (current_frame->frame_refs_short_signaling) {
        // NOTE(zoeliu@google.com):
        //   An example solution for encoder-side implementation on frame refs
        //   short signaling, which is only turned on when the encoder side
        //   decision on ref frames is identical to that at the decoder side.
        current_frame->frame_refs_short_signaling =
            check_frame_refs_short_signaling(cm);
      }

      if (seq_params->order_hint_info.enable_order_hint)
        aom_wb_write_bit(wb, current_frame->frame_refs_short_signaling);

      if (current_frame->frame_refs_short_signaling) {
        const int lst_ref = get_ref_frame_map_idx(cm, LAST_FRAME);
        aom_wb_write_literal(wb, lst_ref, REF_FRAMES_LOG2);

        const int gld_ref = get_ref_frame_map_idx(cm, GOLDEN_FRAME);
        aom_wb_write_literal(wb, gld_ref, REF_FRAMES_LOG2);
      }

      for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
        assert(get_ref_frame_map_idx(cm, ref_frame) != INVALID_IDX);
        if (!current_frame->frame_refs_short_signaling)
          aom_wb_write_literal(wb, get_ref_frame_map_idx(cm, ref_frame),
                               REF_FRAMES_LOG2);
        if (seq_params->frame_id_numbers_present_flag) {
          int i = get_ref_frame_map_idx(cm, ref_frame);
          int frame_id_len = seq_params->frame_id_length;
          int diff_len = seq_params->delta_frame_id_length;
          int delta_frame_id_minus_1 =
              ((cm->current_frame_id - cm->ref_frame_id[i] +
                (1 << frame_id_len)) %
               (1 << frame_id_len)) -
              1;
          if (delta_frame_id_minus_1 < 0 ||
              delta_frame_id_minus_1 >= (1 << diff_len)) {
            aom_internal_error(&cpi->common.error, AOM_CODEC_ERROR,
                               "Invalid delta_frame_id_minus_1");
          }
          aom_wb_write_literal(wb, delta_frame_id_minus_1, diff_len);
        }
      }

      if (!cm->error_resilient_mode && frame_size_override_flag) {
        write_frame_size_with_refs(cm, wb);
      } else {
        write_frame_size(cm, frame_size_override_flag, wb);
      }

      if (!cm->cur_frame_force_integer_mv)
        aom_wb_write_bit(wb, cm->allow_high_precision_mv);
      write_frame_interp_filter(cm->interp_filter, wb);
      aom_wb_write_bit(wb, cm->switchable_motion_mode);
      if (frame_might_allow_ref_frame_mvs(cm)) {
        aom_wb_write_bit(wb, cm->allow_ref_frame_mvs);
      } else {
        assert(cm->allow_ref_frame_mvs == 0);
      }
    }
  }

  const int might_bwd_adapt =
      !(seq_params->reduced_still_picture_hdr) && !(cm->disable_cdf_update);
  if (cm->large_scale_tile)
    assert(cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_DISABLED);

  if (might_bwd_adapt) {
    aom_wb_write_bit(
        wb, cm->refresh_frame_context == REFRESH_FRAME_CONTEXT_DISABLED);
  }

  write_tile_info(cm, saved_wb, wb);
  encode_quantization(cm, wb);
  encode_segmentation(cm, xd, wb);

  const DeltaQInfo *const delta_q_info = &cm->delta_q_info;
  if (delta_q_info->delta_q_present_flag) assert(cm->base_qindex > 0);
  if (cm->base_qindex > 0) {
    aom_wb_write_bit(wb, delta_q_info->delta_q_present_flag);
    if (delta_q_info->delta_q_present_flag) {
      aom_wb_write_literal(wb, get_msb(delta_q_info->delta_q_res), 2);
      xd->current_qindex = cm->base_qindex;
      if (cm->allow_intrabc)
        assert(delta_q_info->delta_lf_present_flag == 0);
      else
        aom_wb_write_bit(wb, delta_q_info->delta_lf_present_flag);
      if (delta_q_info->delta_lf_present_flag) {
        aom_wb_write_literal(wb, get_msb(delta_q_info->delta_lf_res), 2);
        aom_wb_write_bit(wb, delta_q_info->delta_lf_multi);
        av1_reset_loop_filter_delta(xd, av1_num_planes(cm));
      }
    }
  }

  if (cm->all_lossless) {
    assert(!av1_superres_scaled(cm));
  } else {
    if (!cm->coded_lossless) {
      encode_loopfilter(cm, wb);
      encode_cdef(cm, wb);
    }
    encode_restoration_mode(cm, wb);
  }

  // Write TX mode
  if (cm->coded_lossless)
    assert(cm->tx_mode == ONLY_4X4);
  else
    aom_wb_write_bit(wb, cm->tx_mode == TX_MODE_SELECT);

  if (!frame_is_intra_only(cm)) {
    const int use_hybrid_pred =
        current_frame->reference_mode == REFERENCE_MODE_SELECT;

    aom_wb_write_bit(wb, use_hybrid_pred);
  }

  if (current_frame->skip_mode_info.skip_mode_allowed)
    aom_wb_write_bit(wb, current_frame->skip_mode_info.skip_mode_flag);

  if (frame_might_allow_warped_motion(cm))
    aom_wb_write_bit(wb, cm->allow_warped_motion);
  else
    assert(!cm->allow_warped_motion);

  aom_wb_write_bit(wb, cm->reduced_tx_set_used);

  if (!frame_is_intra_only(cm)) write_global_motion(cpi, wb);

  if (seq_params->film_grain_params_present &&
      (cm->show_frame || cm->showable_frame))
    write_film_grain_params(cpi, wb);

  if (cm->large_scale_tile) write_ext_tile_info(cm, saved_wb, wb);
}

static int choose_size_bytes(uint32_t size, int spare_msbs) {
  // Choose the number of bytes required to represent size, without
  // using the 'spare_msbs' number of most significant bits.

  // Make sure we will fit in 4 bytes to start with..
  if (spare_msbs > 0 && size >> (32 - spare_msbs) != 0) return -1;

  // Normalise to 32 bits
  size <<= spare_msbs;

  if (size >> 24 != 0)
    return 4;
  else if (size >> 16 != 0)
    return 3;
  else if (size >> 8 != 0)
    return 2;
  else
    return 1;
}

static void mem_put_varsize(uint8_t *const dst, const int sz, const int val) {
  switch (sz) {
    case 1: dst[0] = (uint8_t)(val & 0xff); break;
    case 2: mem_put_le16(dst, val); break;
    case 3: mem_put_le24(dst, val); break;
    case 4: mem_put_le32(dst, val); break;
    default: assert(0 && "Invalid size"); break;
  }
}

static int remux_tiles(const AV1_COMMON *const cm, uint8_t *dst,
                       const uint32_t data_size, const uint32_t max_tile_size,
                       const uint32_t max_tile_col_size,
                       int *const tile_size_bytes,
                       int *const tile_col_size_bytes) {
  // Choose the tile size bytes (tsb) and tile column size bytes (tcsb)
  int tsb;
  int tcsb;

  if (cm->large_scale_tile) {
    // The top bit in the tile size field indicates tile copy mode, so we
    // have 1 less bit to code the tile size
    tsb = choose_size_bytes(max_tile_size, 1);
    tcsb = choose_size_bytes(max_tile_col_size, 0);
  } else {
    tsb = choose_size_bytes(max_tile_size, 0);
    tcsb = 4;  // This is ignored
    (void)max_tile_col_size;
  }

  assert(tsb > 0);
  assert(tcsb > 0);

  *tile_size_bytes = tsb;
  *tile_col_size_bytes = tcsb;
  if (tsb == 4 && tcsb == 4) return data_size;

  uint32_t wpos = 0;
  uint32_t rpos = 0;

  if (cm->large_scale_tile) {
    int tile_row;
    int tile_col;

    for (tile_col = 0; tile_col < cm->tile_cols; tile_col++) {
      // All but the last column has a column header
      if (tile_col < cm->tile_cols - 1) {
        uint32_t tile_col_size = mem_get_le32(dst + rpos);
        rpos += 4;

        // Adjust the tile column size by the number of bytes removed
        // from the tile size fields.
        tile_col_size -= (4 - tsb) * cm->tile_rows;

        mem_put_varsize(dst + wpos, tcsb, tile_col_size);
        wpos += tcsb;
      }

      for (tile_row = 0; tile_row < cm->tile_rows; tile_row++) {
        // All, including the last row has a header
        uint32_t tile_header = mem_get_le32(dst + rpos);
        rpos += 4;

        // If this is a copy tile, we need to shift the MSB to the
        // top bit of the new width, and there is no data to copy.
        if (tile_header >> 31 != 0) {
          if (tsb < 4) tile_header >>= 32 - 8 * tsb;
          mem_put_varsize(dst + wpos, tsb, tile_header);
          wpos += tsb;
        } else {
          mem_put_varsize(dst + wpos, tsb, tile_header);
          wpos += tsb;

          tile_header += AV1_MIN_TILE_SIZE_BYTES;
          memmove(dst + wpos, dst + rpos, tile_header);
          rpos += tile_header;
          wpos += tile_header;
        }
      }
    }

    assert(rpos > wpos);
    assert(rpos == data_size);

    return wpos;
  }
  const int n_tiles = cm->tile_cols * cm->tile_rows;
  int n;

  for (n = 0; n < n_tiles; n++) {
    int tile_size;

    if (n == n_tiles - 1) {
      tile_size = data_size - rpos;
    } else {
      tile_size = mem_get_le32(dst + rpos);
      rpos += 4;
      mem_put_varsize(dst + wpos, tsb, tile_size);
      tile_size += AV1_MIN_TILE_SIZE_BYTES;
      wpos += tsb;
    }

    memmove(dst + wpos, dst + rpos, tile_size);

    rpos += tile_size;
    wpos += tile_size;
  }

  assert(rpos > wpos);
  assert(rpos == data_size);

  return wpos;
}

uint32_t av1_write_obu_header(AV1_COMP *const cpi, OBU_TYPE obu_type,
                              int obu_extension, uint8_t *const dst) {
  if (cpi->keep_level_stats &&
      (obu_type == OBU_FRAME || obu_type == OBU_FRAME_HEADER))
    ++cpi->frame_header_count;

  struct aom_write_bit_buffer wb = { dst, 0 };
  uint32_t size = 0;

  aom_wb_write_literal(&wb, 0, 1);  // forbidden bit.
  aom_wb_write_literal(&wb, (int)obu_type, 4);
  aom_wb_write_literal(&wb, obu_extension ? 1 : 0, 1);
  aom_wb_write_literal(&wb, 1, 1);  // obu_has_payload_length_field
  aom_wb_write_literal(&wb, 0, 1);  // reserved

  if (obu_extension) {
    aom_wb_write_literal(&wb, obu_extension & 0xFF, 8);
  }

  size = aom_wb_bytes_written(&wb);
  return size;
}

int av1_write_uleb_obu_size(uint32_t obu_header_size, uint32_t obu_payload_size,
                            uint8_t *dest) {
  const uint32_t obu_size = obu_payload_size;
  const uint32_t offset = obu_header_size;
  size_t coded_obu_size = 0;

  if (aom_uleb_encode(obu_size, sizeof(obu_size), dest + offset,
                      &coded_obu_size) != 0) {
    return AOM_CODEC_ERROR;
  }

  return AOM_CODEC_OK;
}

static size_t obu_memmove(uint32_t obu_header_size, uint32_t obu_payload_size,
                          uint8_t *data) {
  const size_t length_field_size = aom_uleb_size_in_bytes(obu_payload_size);
  const uint32_t move_dst_offset =
      (uint32_t)length_field_size + obu_header_size;
  const uint32_t move_src_offset = obu_header_size;
  const uint32_t move_size = obu_payload_size;
  memmove(data + move_dst_offset, data + move_src_offset, move_size);
  return length_field_size;
}

static void add_trailing_bits(struct aom_write_bit_buffer *wb) {
  if (aom_wb_is_byte_aligned(wb)) {
    aom_wb_write_literal(wb, 0x80, 8);
  } else {
    // assumes that the other bits are already 0s
    aom_wb_write_bit(wb, 1);
  }
}

static void write_bitstream_level(AV1_LEVEL seq_level_idx,
                                  struct aom_write_bit_buffer *wb) {
  assert(is_valid_seq_level_idx(seq_level_idx));
  aom_wb_write_literal(wb, seq_level_idx, LEVEL_BITS);
}

uint32_t av1_write_sequence_header_obu(AV1_COMP *cpi, uint8_t *const dst) {
  AV1_COMMON *const cm = &cpi->common;
  struct aom_write_bit_buffer wb = { dst, 0 };
  uint32_t size = 0;

  write_profile(cm->seq_params.profile, &wb);

  // Still picture or not
  aom_wb_write_bit(&wb, cm->seq_params.still_picture);
  assert(IMPLIES(!cm->seq_params.still_picture,
                 !cm->seq_params.reduced_still_picture_hdr));
  // whether to use reduced still picture header
  aom_wb_write_bit(&wb, cm->seq_params.reduced_still_picture_hdr);

  if (cm->seq_params.reduced_still_picture_hdr) {
    assert(cm->timing_info_present == 0);
    assert(cm->seq_params.decoder_model_info_present_flag == 0);
    assert(cm->seq_params.display_model_info_present_flag == 0);
    write_bitstream_level(cm->seq_params.seq_level_idx[0], &wb);
  } else {
    aom_wb_write_bit(&wb, cm->timing_info_present);  // timing info present flag

    if (cm->timing_info_present) {
      // timing_info
      write_timing_info_header(cm, &wb);
      aom_wb_write_bit(&wb, cm->seq_params.decoder_model_info_present_flag);
      if (cm->seq_params.decoder_model_info_present_flag) {
        write_decoder_model_info(cm, &wb);
      }
    }
    aom_wb_write_bit(&wb, cm->seq_params.display_model_info_present_flag);
    aom_wb_write_literal(&wb, cm->seq_params.operating_points_cnt_minus_1,
                         OP_POINTS_CNT_MINUS_1_BITS);
    int i;
    for (i = 0; i < cm->seq_params.operating_points_cnt_minus_1 + 1; i++) {
      aom_wb_write_literal(&wb, cm->seq_params.operating_point_idc[i],
                           OP_POINTS_IDC_BITS);
      write_bitstream_level(cm->seq_params.seq_level_idx[i], &wb);
      if (cm->seq_params.seq_level_idx[i] >= SEQ_LEVEL_4_0)
        aom_wb_write_bit(&wb, cm->seq_params.tier[i]);
      if (cm->seq_params.decoder_model_info_present_flag) {
        aom_wb_write_bit(&wb,
                         cm->op_params[i].decoder_model_param_present_flag);
        if (cm->op_params[i].decoder_model_param_present_flag)
          write_dec_model_op_parameters(cm, &wb, i);
      }
      if (cm->seq_params.display_model_info_present_flag) {
        aom_wb_write_bit(&wb,
                         cm->op_params[i].display_model_param_present_flag);
        if (cm->op_params[i].display_model_param_present_flag) {
          assert(cm->op_params[i].initial_display_delay <= 10);
          aom_wb_write_literal(&wb, cm->op_params[i].initial_display_delay - 1,
                               4);
        }
      }
    }
  }
  write_sequence_header(&cm->seq_params, &wb);

  write_color_config(&cm->seq_params, &wb);

  aom_wb_write_bit(&wb, cm->seq_params.film_grain_params_present);

  add_trailing_bits(&wb);

  size = aom_wb_bytes_written(&wb);
  return size;
}

static uint32_t write_frame_header_obu(AV1_COMP *cpi,
                                       struct aom_write_bit_buffer *saved_wb,
                                       uint8_t *const dst,
                                       int append_trailing_bits) {
  struct aom_write_bit_buffer wb = { dst, 0 };
  write_uncompressed_header_obu(cpi, saved_wb, &wb);
  if (append_trailing_bits) add_trailing_bits(&wb);
  return aom_wb_bytes_written(&wb);
}

static uint32_t write_tile_group_header(uint8_t *const dst, int start_tile,
                                        int end_tile, int tiles_log2,
                                        int tile_start_and_end_present_flag) {
  struct aom_write_bit_buffer wb = { dst, 0 };
  uint32_t size = 0;

  if (!tiles_log2) return size;

  aom_wb_write_bit(&wb, tile_start_and_end_present_flag);

  if (tile_start_and_end_present_flag) {
    aom_wb_write_literal(&wb, start_tile, tiles_log2);
    aom_wb_write_literal(&wb, end_tile, tiles_log2);
  }

  size = aom_wb_bytes_written(&wb);
  return size;
}

typedef struct {
  uint8_t *frame_header;
  size_t obu_header_byte_offset;
  size_t total_length;
} FrameHeaderInfo;

static uint32_t write_tiles_in_tg_obus(AV1_COMP *const cpi, uint8_t *const dst,
                                       struct aom_write_bit_buffer *saved_wb,
                                       uint8_t obu_extension_header,
                                       const FrameHeaderInfo *fh_info,
                                       int *const largest_tile_id) {
  AV1_COMMON *const cm = &cpi->common;
  aom_writer mode_bc;
  int tile_row, tile_col;
  // Store the location and size of each tile's data in the bitstream:
  TileBufferEnc tile_buffers[MAX_TILE_ROWS][MAX_TILE_COLS];
  uint32_t total_size = 0;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;
  unsigned int tile_size = 0;
  unsigned int max_tile_size = 0;
  unsigned int max_tile_col_size = 0;
  const int n_log2_tiles = cm->log2_tile_rows + cm->log2_tile_cols;
  // Fixed size tile groups for the moment
  const int num_tg_hdrs = cm->num_tg;
  const int tg_size =
      (cm->large_scale_tile)
          ? 1
          : (tile_rows * tile_cols + num_tg_hdrs - 1) / num_tg_hdrs;
  int tile_count = 0;
  int curr_tg_data_size = 0;
  uint8_t *data = dst;
  int new_tg = 1;
  const int have_tiles = tile_cols * tile_rows > 1;
  int first_tg = 1;

  *largest_tile_id = 0;

  if (cm->large_scale_tile) {
    // For large_scale_tile case, we always have only one tile group, so it can
    // be written as an OBU_FRAME.
    const OBU_TYPE obu_type = OBU_FRAME;
    const uint32_t tg_hdr_size = av1_write_obu_header(cpi, obu_type, 0, data);
    data += tg_hdr_size;

    const uint32_t frame_header_size =
        write_frame_header_obu(cpi, saved_wb, data, 0);
    data += frame_header_size;
    total_size += frame_header_size;

#define EXT_TILE_DEBUG 0
#if EXT_TILE_DEBUG
    {
      char fn[20] = "./fh";
      fn[4] = cm->current_frame.frame_number / 100 + '0';
      fn[5] = (cm->current_frame.frame_number % 100) / 10 + '0';
      fn[6] = (cm->current_frame.frame_number % 10) + '0';
      fn[7] = '\0';
      av1_print_uncompressed_frame_header(data - frame_header_size,
                                          frame_header_size, fn);
    }
#endif  // EXT_TILE_DEBUG
#undef EXT_TILE_DEBUG

    int tile_size_bytes = 0;
    int tile_col_size_bytes = 0;

    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      TileInfo tile_info;
      const int is_last_col = (tile_col == tile_cols - 1);
      const uint32_t col_offset = total_size;

      av1_tile_set_col(&tile_info, cm, tile_col);

      // The last column does not have a column header
      if (!is_last_col) total_size += 4;

      for (tile_row = 0; tile_row < tile_rows; tile_row++) {
        TileBufferEnc *const buf = &tile_buffers[tile_row][tile_col];
        const int data_offset = have_tiles ? 4 : 0;
        const int tile_idx = tile_row * tile_cols + tile_col;
        TileDataEnc *this_tile = &cpi->tile_data[tile_idx];
        av1_tile_set_row(&tile_info, cm, tile_row);

        buf->data = dst + total_size + tg_hdr_size;

        // Is CONFIG_EXT_TILE = 1, every tile in the row has a header,
        // even for the last one, unless no tiling is used at all.
        total_size += data_offset;
        cpi->td.mb.e_mbd.tile_ctx = &this_tile->tctx;
        mode_bc.allow_update_cdf = !cm->large_scale_tile;
        mode_bc.allow_update_cdf =
            mode_bc.allow_update_cdf && !cm->disable_cdf_update;
        aom_start_encode(&mode_bc, buf->data + data_offset);
        write_modes(cpi, &tile_info, &mode_bc, tile_row, tile_col);
        aom_stop_encode(&mode_bc);
        tile_size = mode_bc.pos;
        buf->size = tile_size;

        // Record the maximum tile size we see, so we can compact headers later.
        if (tile_size > max_tile_size) {
          max_tile_size = tile_size;
          *largest_tile_id = tile_cols * tile_row + tile_col;
        }

        if (have_tiles) {
          // tile header: size of this tile, or copy offset
          uint32_t tile_header = tile_size - AV1_MIN_TILE_SIZE_BYTES;
          const int tile_copy_mode =
              ((AOMMAX(cm->tile_width, cm->tile_height) << MI_SIZE_LOG2) <= 256)
                  ? 1
                  : 0;

          // If tile_copy_mode = 1, check if this tile is a copy tile.
          // Very low chances to have copy tiles on the key frames, so don't
          // search on key frames to reduce unnecessary search.
          if (cm->current_frame.frame_type != KEY_FRAME && tile_copy_mode) {
            const int identical_tile_offset =
                find_identical_tile(tile_row, tile_col, tile_buffers);

            // Indicate a copy-tile by setting the most significant bit.
            // The row-offset to copy from is stored in the highest byte.
            // remux_tiles will move these around later
            if (identical_tile_offset > 0) {
              tile_size = 0;
              tile_header = identical_tile_offset | 0x80;
              tile_header <<= 24;
            }
          }

          mem_put_le32(buf->data, tile_header);
        }

        total_size += tile_size;
      }

      if (!is_last_col) {
        uint32_t col_size = total_size - col_offset - 4;
        mem_put_le32(dst + col_offset + tg_hdr_size, col_size);

        // Record the maximum tile column size we see.
        max_tile_col_size = AOMMAX(max_tile_col_size, col_size);
      }
    }

    if (have_tiles) {
      total_size = remux_tiles(cm, data, total_size - frame_header_size,
                               max_tile_size, max_tile_col_size,
                               &tile_size_bytes, &tile_col_size_bytes);
      total_size += frame_header_size;
    }

    // In EXT_TILE case, only use 1 tile group. Follow the obu syntax, write
    // current tile group size before tile data(include tile column header).
    // Tile group size doesn't include the bytes storing tg size.
    total_size += tg_hdr_size;
    const uint32_t obu_payload_size = total_size - tg_hdr_size;
    const size_t length_field_size =
        obu_memmove(tg_hdr_size, obu_payload_size, dst);
    if (av1_write_uleb_obu_size(tg_hdr_size, obu_payload_size, dst) !=
        AOM_CODEC_OK) {
      assert(0);
    }
    total_size += (uint32_t)length_field_size;
    saved_wb->bit_buffer += length_field_size;

    // Now fill in the gaps in the uncompressed header.
    if (have_tiles) {
      assert(tile_col_size_bytes >= 1 && tile_col_size_bytes <= 4);
      aom_wb_overwrite_literal(saved_wb, tile_col_size_bytes - 1, 2);

      assert(tile_size_bytes >= 1 && tile_size_bytes <= 4);
      aom_wb_overwrite_literal(saved_wb, tile_size_bytes - 1, 2);
    }
    return total_size;
  }

  uint32_t obu_header_size = 0;
  uint8_t *tile_data_start = dst + total_size;
  for (tile_row = 0; tile_row < tile_rows; tile_row++) {
    TileInfo tile_info;
    av1_tile_set_row(&tile_info, cm, tile_row);

    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      const int tile_idx = tile_row * tile_cols + tile_col;
      TileBufferEnc *const buf = &tile_buffers[tile_row][tile_col];
      TileDataEnc *this_tile = &cpi->tile_data[tile_idx];
      int is_last_tile_in_tg = 0;

      if (new_tg) {
        data = dst + total_size;

        // A new tile group begins at this tile.  Write the obu header and
        // tile group header
        const OBU_TYPE obu_type =
            (num_tg_hdrs == 1) ? OBU_FRAME : OBU_TILE_GROUP;
        curr_tg_data_size =
            av1_write_obu_header(cpi, obu_type, obu_extension_header, data);
        obu_header_size = curr_tg_data_size;

        if (num_tg_hdrs == 1) {
          curr_tg_data_size += write_frame_header_obu(
              cpi, saved_wb, data + curr_tg_data_size, 0);
        }
        curr_tg_data_size += write_tile_group_header(
            data + curr_tg_data_size, tile_idx,
            AOMMIN(tile_idx + tg_size - 1, tile_cols * tile_rows - 1),
            n_log2_tiles, cm->num_tg > 1);
        total_size += curr_tg_data_size;
        tile_data_start += curr_tg_data_size;
        new_tg = 0;
        tile_count = 0;
      }
      tile_count++;
      av1_tile_set_col(&tile_info, cm, tile_col);

      if (tile_count == tg_size || tile_idx == (tile_cols * tile_rows - 1)) {
        is_last_tile_in_tg = 1;
        new_tg = 1;
      } else {
        is_last_tile_in_tg = 0;
      }

      buf->data = dst + total_size;

      // The last tile of the tile group does not have a header.
      if (!is_last_tile_in_tg) total_size += 4;

      cpi->td.mb.e_mbd.tile_ctx = &this_tile->tctx;
      mode_bc.allow_update_cdf = 1;
      mode_bc.allow_update_cdf =
          mode_bc.allow_update_cdf && !cm->disable_cdf_update;
      const int num_planes = av1_num_planes(cm);
      av1_reset_loop_restoration(&cpi->td.mb.e_mbd, num_planes);

      aom_start_encode(&mode_bc, dst + total_size);
      write_modes(cpi, &tile_info, &mode_bc, tile_row, tile_col);
      aom_stop_encode(&mode_bc);
      tile_size = mode_bc.pos;
      assert(tile_size >= AV1_MIN_TILE_SIZE_BYTES);

      curr_tg_data_size += (tile_size + (is_last_tile_in_tg ? 0 : 4));
      buf->size = tile_size;
      if (tile_size > max_tile_size) {
        *largest_tile_id = tile_cols * tile_row + tile_col;
        max_tile_size = tile_size;
      }

      if (!is_last_tile_in_tg) {
        // size of this tile
        mem_put_le32(buf->data, tile_size - AV1_MIN_TILE_SIZE_BYTES);
      } else {
        // write current tile group size
        const uint32_t obu_payload_size = curr_tg_data_size - obu_header_size;
        const size_t length_field_size =
            obu_memmove(obu_header_size, obu_payload_size, data);
        if (av1_write_uleb_obu_size(obu_header_size, obu_payload_size, data) !=
            AOM_CODEC_OK) {
          assert(0);
        }
        curr_tg_data_size += (int)length_field_size;
        total_size += (uint32_t)length_field_size;
        tile_data_start += length_field_size;
        if (num_tg_hdrs == 1) {
          // if this tg is combined with the frame header then update saved
          // frame header base offset accroding to length field size
          saved_wb->bit_buffer += length_field_size;
        }

        if (!first_tg && cm->error_resilient_mode) {
          // Make room for a duplicate Frame Header OBU.
          memmove(data + fh_info->total_length, data, curr_tg_data_size);

          // Insert a copy of the Frame Header OBU.
          memcpy(data, fh_info->frame_header, fh_info->total_length);

          // Force context update tile to be the first tile in error
          // resiliant mode as the duplicate frame headers will have
          // context_update_tile_id set to 0
          *largest_tile_id = 0;

          // Rewrite the OBU header to change the OBU type to Redundant Frame
          // Header.
          av1_write_obu_header(cpi, OBU_REDUNDANT_FRAME_HEADER,
                               obu_extension_header,
                               &data[fh_info->obu_header_byte_offset]);

          data += fh_info->total_length;

          curr_tg_data_size += (int)(fh_info->total_length);
          total_size += (uint32_t)(fh_info->total_length);
        }
        first_tg = 0;
      }

      total_size += tile_size;
    }
  }

  if (have_tiles) {
    // Fill in context_update_tile_id indicating the tile to use for the
    // cdf update. The encoder currently sets it to the largest tile
    // (but is up to the encoder)
    aom_wb_overwrite_literal(saved_wb, *largest_tile_id,
                             cm->log2_tile_cols + cm->log2_tile_rows);
    // If more than one tile group. tile_size_bytes takes the default value 4
    // and does not need to be set. For a single tile group it is set in the
    // section below.
    if (num_tg_hdrs == 1) {
      int tile_size_bytes = 4, unused;
      const uint32_t tile_data_offset = (uint32_t)(tile_data_start - dst);
      const uint32_t tile_data_size = total_size - tile_data_offset;

      total_size =
          remux_tiles(cm, tile_data_start, tile_data_size, max_tile_size,
                      max_tile_col_size, &tile_size_bytes, &unused);
      total_size += tile_data_offset;
      assert(tile_size_bytes >= 1 && tile_size_bytes <= 4);

      aom_wb_overwrite_literal(saved_wb, tile_size_bytes - 1, 2);

      // Update the OBU length if remux_tiles() reduced the size.
      uint64_t payload_size;
      size_t length_field_size;
      int res =
          aom_uleb_decode(dst + obu_header_size, total_size - obu_header_size,
                          &payload_size, &length_field_size);
      assert(res == 0);
      (void)res;

      const uint64_t new_payload_size =
          total_size - obu_header_size - length_field_size;
      if (new_payload_size != payload_size) {
        size_t new_length_field_size;
        res = aom_uleb_encode(new_payload_size, length_field_size,
                              dst + obu_header_size, &new_length_field_size);
        assert(res == 0);
        if (new_length_field_size < length_field_size) {
          const size_t src_offset = obu_header_size + length_field_size;
          const size_t dst_offset = obu_header_size + new_length_field_size;
          memmove(dst + dst_offset, dst + src_offset, (size_t)payload_size);
          total_size -= (int)(length_field_size - new_length_field_size);
        }
      }
    }
  }
  return total_size;
}

int av1_pack_bitstream(AV1_COMP *const cpi, uint8_t *dst, size_t *size,
                       int *const largest_tile_id) {
  uint8_t *data = dst;
  uint32_t data_size;
  AV1_COMMON *const cm = &cpi->common;
  uint32_t obu_header_size = 0;
  uint32_t obu_payload_size = 0;
  FrameHeaderInfo fh_info = { NULL, 0, 0 };
  const uint8_t obu_extension_header =
      cm->temporal_layer_id << 5 | cm->spatial_layer_id << 3 | 0;

  // If no non-zero delta_q has been used, reset delta_q_present_flag
  if (cm->delta_q_info.delta_q_present_flag && cpi->deltaq_used == 0) {
    cm->delta_q_info.delta_q_present_flag = 0;
  }

#if CONFIG_BITSTREAM_DEBUG
  bitstream_queue_reset_write();
#endif

  cpi->frame_header_count = 0;

  // The TD is now written outside the frame encode loop

  // write sequence header obu if KEY_FRAME, preceded by 4-byte size
  if (cm->current_frame.frame_type == KEY_FRAME && cm->show_frame) {
    obu_header_size = av1_write_obu_header(cpi, OBU_SEQUENCE_HEADER, 0, data);

    obu_payload_size =
        av1_write_sequence_header_obu(cpi, data + obu_header_size);
    const size_t length_field_size =
        obu_memmove(obu_header_size, obu_payload_size, data);
    if (av1_write_uleb_obu_size(obu_header_size, obu_payload_size, data) !=
        AOM_CODEC_OK) {
      return AOM_CODEC_ERROR;
    }

    data += obu_header_size + obu_payload_size + length_field_size;
  }

  const int write_frame_header =
      (cm->num_tg > 1 || encode_show_existing_frame(cm));
  struct aom_write_bit_buffer saved_wb;
  if (write_frame_header) {
    // Write Frame Header OBU.
    fh_info.frame_header = data;
    obu_header_size =
        av1_write_obu_header(cpi, OBU_FRAME_HEADER, obu_extension_header, data);
    obu_payload_size =
        write_frame_header_obu(cpi, &saved_wb, data + obu_header_size, 1);

    const size_t length_field_size =
        obu_memmove(obu_header_size, obu_payload_size, data);
    if (av1_write_uleb_obu_size(obu_header_size, obu_payload_size, data) !=
        AOM_CODEC_OK) {
      return AOM_CODEC_ERROR;
    }

    fh_info.obu_header_byte_offset = 0;
    fh_info.total_length =
        obu_header_size + obu_payload_size + length_field_size;
    data += fh_info.total_length;

    // Since length_field_size is determined adaptively after frame header
    // encoding, saved_wb must be adjusted accordingly.
    saved_wb.bit_buffer += length_field_size;
  }

  if (encode_show_existing_frame(cm)) {
    data_size = 0;
  } else {
    //  Each tile group obu will be preceded by 4-byte size of the tile group
    //  obu
    data_size = write_tiles_in_tg_obus(
        cpi, data, &saved_wb, obu_extension_header, &fh_info, largest_tile_id);
  }
  data += data_size;
  *size = data - dst;
  return AOM_CODEC_OK;
}
