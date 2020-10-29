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

#ifndef AOM_AV1_ENCODER_AV1_ML_PARTITION_MODELS_H_
#define AOM_AV1_ENCODER_AV1_ML_PARTITION_MODELS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/ml.h"

// TODO(kyslov): Replace with proper weights after training AV1 models

#define FEATURES 4
/*
static const float av1_var_part_nn_weights_64_layer0[FEATURES * 8] = {
  0.35755366f,  0.86281112f,  -0.20871686f, 0.0409634f,   0.97305766f,
  0.75510254f,  0.04860447f,  0.77095283f,  -0.44105278f, -0.3755049f,
  -0.08456618f, 1.1821136f,   -0.73956301f, 1.30016453f,  0.45566902f,
  0.4742967f,   0.44213975f,  0.4876028f,   0.26720522f,  -0.34429858f,
  -0.25148252f, -0.49623932f, -0.46747941f, -0.36656624f, 0.10213375f,
  0.60262819f,  -0.54788715f, -0.27272022f, 1.0995462f,   -0.36338376f,
  -0.64836313f, 0.16057039f,  1.02782791f,  0.9985311f,   0.90607883f,
  0.80570411f,  -0.07750863f, -0.74006402f, 1.72839526f,  1.72355343f,
  1.69288916f,  1.59102043f,  0.14140216f,  -1.47262839f, 0.4262519f,
  -0.33805936f, -0.02449707f, 0.67203692f
};

static const float av1_var_part_nn_bias_64_layer0[8] = {
  0.39995694f, 0.65593756f, 1.12876737f,  1.28790576f,
  0.53468556f, 0.3177908f,  -0.74388266f, -1.81131248f
};

static const float av1_var_part_nn_weights_64_layer1[8] = {
  -1.31174053f, 0.69696917f, 0.78721456f, 0.45326379f,
  0.79258322f,  1.74626188f, -5.41831f,   3.33887435f
};

static const float av1_var_part_nn_bias_64_layer1[1] = { -0.90951047f };

static const float av1_var_part_means_64[FEATURES] = {
  5.36750249f, 11.58023127f, 0.25550964f, 0.23809917f, 0.24650665f, 0.22117687f
};
static const float av1_var_part_vars_64[FEATURES] = {
  0.89599769f, 2.2686018f, 0.02568608f, 0.02523411f, 0.02443085f, 0.01922085f
};
*/

/* 6-features full ------
static const float av1_var_part_nn_weights_64_layer0[FEATURES * 8] = {
  -8.54925369e-01f, 7.97433626e-01f, -3.73334232e-01f, -3.54351323e-01f,
-4.17120892e-01f, -3.10031378e-01f, 7.77357747e-01f,
-1.55523059e+00f, 3.80701829e-01f, 5.18489285e-01f, 2.64749615e-01f,
-9.91847162e-02, -3.31459929e-02f, -1.82533077e-01f,
3.26929899e-01f, 3.77273646e-01f, -5.66529279e-01f, 1.09873552e+00f,
  -2.23829340e-04f, -3.64472628e-01f, 1.24501263e+00f,
1.28083112e+00f, 1.36294960e+00f, 1.31340377e+00f,
  -9.92732584e-01f, 9.78253815e-01f, -1.90091364e-01f, -2.58895923e-01f,
-1.03148528e-01f, -3.82784003e-01f, 2.65331219e-01f,
-9.69803362e-01f, 4.08802385e-01f, 7.23463259e-01f,
1.07728872e+00f, 1.11206991e+00f, 1.80899210e-01f, -8.34760680e-01f,
-1.04057362e+00f, -8.54278068e-02f, 1.08678103e-01f, 1.82654318e-01f,
  5.45675703e-01f, -1.00136355e+00f, 1.08327852e+00f,
1.14782434e+00f, 1.11757351e+00f, 1.27903610e+00f
};

static const float av1_var_part_nn_bias_64_layer0[8] = {
  -0.37609937f, -0.60197745f, -0.20067775f, -0.38069221f,
  0.6993178f, -0.94040267f, -1.43647201f, -0.02289873f
};

static const float av1_var_part_nn_weights_64_layer1[8] = {
  -1.4932489f, 1.32700348f, 0.91740533f, -4.43162268f, 1.49120542f, 3.31915758f,
-1.98026192f, -2.88562596f
};

static const float av1_var_part_nn_bias_64_layer1[1] = { 0.40843818 };

static const float av1_var_part_means_64[FEATURES] = {
  6.28225573f, 12.15998333f,  0.24116005f,  0.23321305f, 0.25283904f, 0.2421259f
};
static const float av1_var_part_vars_64[FEATURES] = {
  1.30489428f, 2.87956562f, 0.02229749f, 0.02347817f, 0.02518943f, 0.02436969f
};
*/
/* 4 features - full =======
static const float av1_var_part_nn_weights_64_layer0[FEATURES * 8] = {
  1.49263025f, -1.04689216f, 0.11668445f, 0.3774183f,
  -1.33524408f, -0.94913651f, 0.28547942f, -0.25631165f,
  -0.53522048f, -0.0678052f, -1.80100425f, 0.19909253f,
  -0.60001088f, 1.86405481f, -0.11653572f, 0.55590594f,
  -1.63683194f, 2.9963502f, -0.06310161f, 0.64628252f,
  0.49737501f, -1.34459763f, -1.49205269f, -0.75898519f,
  -1.01670898f, 0.63320291f, 0.40182596f, 1.78388687f,
  -1.78042721f, 2.17986441f, 0.05345279f, 0.15334089f
};

static const float av1_var_part_nn_bias_64_layer0[8] = {
  0.83759125f, -1.78475529f,  1.48752187f,  1.94559704f,
-0.47844072f, 1.45196338f,  1.36688362f,  1.05347522f
};

static const float av1_var_part_nn_weights_64_layer1[8] = {
  -0.41295042f, -0.53543211f, -0.65871394f, 0.91720965f, -1.0951356f,
0.78317512f, 0.52407166f, 1.07183618f
};

static const float av1_var_part_nn_bias_64_layer1[1] = { -1.51922838f };

static const float av1_var_part_means_64[FEATURES] = {
  6.12349212f, 12.05651136f, 0.10803822f, 0.43453482f
};
static const float av1_var_part_vars_64[FEATURES] = {
  1.3949857f, 3.0868353f, 0.00416123f, 0.01896546f
};
*/
/* 4 features full - rdo
static const float av1_var_part_nn_weights_64_layer0[FEATURES * 8] = {
  -1.27279737f, 2.29205384f, 0.16531761f, 2.06375128f,
  -0.77134797f, 1.13970645f, -1.46291045f, 0.74321657f,
  -2.38735165f, 2.90485219f, 0.27194302f, 0.35763905f,
  -0.71647142f, 0.03557003f, 0.32614374f, -0.74569224f,
  -1.13890672f, 2.45435897f, -1.27462625f, -1.0449622f,
  -1.14151541f, -0.14639661f, 0.01106748f, 1.29590258f,
  -1.60081949f, 1.89189825f, -2.79840838f, 1.21301515f,
  -0.27875942f, 1.96411103f, 1.66237571f, 0.40433006f
};

static const float av1_var_part_nn_bias_64_layer0[8] = {
  2.58547305f, 2.44019685f, 1.31949208f, -0.43862397f, 2.53517553f, 0.34355691f,
-4.08708222f,  1.0794201f
};

static const float av1_var_part_nn_weights_64_layer1[8] = {
  0.64468436f, -0.65440742f, 0.36298593f, -0.20791952f, 0.50178626f,
0.32535124f, -0.55393106f, -0.38672159f
};

static const float av1_var_part_nn_bias_64_layer1[1] = { -0.26285466f };

static const float av1_var_part_means_64[FEATURES] = {
  6.20310378f, 12.0040552f, 0.10718053f, 0.43750866f
};

static const float av1_var_part_vars_64[FEATURES] = {
  1.43159293f, 3.12349129f, 0.00428601f, 0.01953825f
};
*/

static const float av1_var_part_nn_weights_64_layer0[FEATURES * 8] = {
  0.96458016f,  -2.16374735f, 0.43554644f,  0.27604013f,  -1.17570017f,
  0.85239039f,  -1.847042f,   0.82519525f,  -1.95331647f, 2.46487279f,
  -0.16775882f, 0.33033214f,  0.04853412f,  -0.10376066f, -1.32769625f,
  -0.09650413f, 0.52757402f,  -2.17009251f, -0.50047693f, -1.99055089f,
  -1.8508833f,  0.05680268f,  0.67704912f,  1.12248542f,  0.48833869f,
  0.41930461f,  -0.29586123f, 1.07133729f,  0.60099965f,  1.43858737f,
  0.00477526f,  1.11522038f
};

static const float av1_var_part_nn_bias_64_layer0[8] = {
  2.14801393f,  -3.50417954f, 1.5887147f,  0.91487323f,
  -3.74198426f, -0.46068389f, 2.98338532f, 1.51120471f
};

static const float av1_var_part_nn_weights_64_layer1[8] = {
  -0.86394967f, -0.54303225f, 0.62828063f, -0.40975028f,
  1.58661187f,  0.35147492f,  0.82013612f, -0.56718946f
};

static const float av1_var_part_nn_bias_64_layer1[1] = { 0.06911727f };

static const float av1_var_part_means_64[FEATURES] = {
  5.92435081f, 11.82305163f, 0.11360111f, 0.42826822f
};

static const float av1_var_part_vars_64[FEATURES] = { 1.62618851f, 2.51088417f,
                                                      0.00399878f,
                                                      0.01717729f };

static const NN_CONFIG av1_var_part_nnconfig_64 = {
  FEATURES,  // num_inputs
  1,         // num_outputs
  1,         // num_hidden_layers
  {
      8,
  },  // num_hidden_nodes
  {
      av1_var_part_nn_weights_64_layer0,
      av1_var_part_nn_weights_64_layer1,
  },
  {
      av1_var_part_nn_bias_64_layer0,
      av1_var_part_nn_bias_64_layer1,
  },
};

/*
static const float av1_var_part_nn_weights_32_layer0[FEATURES * 8] = {
  0.97886049f,  -1.66262011f, 0.94902798f,  0.7080922f,   0.91181186f,
  0.35222601f,  -0.04428585f, 0.42086472f,  -0.0206325f,  -0.77937809f,
  -0.70947522f, -1.24463119f, 0.23739497f,  -1.34327359f, 0.01024804f,
  0.4544633f,   -0.96907661f, 0.67279522f,  0.23180693f,  1.54063368f,
  -0.15700707f, 0.18597331f,  0.34167589f,  0.40736558f,  0.69213366f,
  -1.33584593f, 1.21190814f,  1.26725267f,  1.21284802f,  1.26611399f,
  0.17546514f,  -0.30248399f, -1.32589316f, -1.37432674f, -1.37423023f,
  -1.26890855f, 0.12166347f,  -0.94565678f, -1.47475267f, -0.69279948f,
  -0.10166587f, -0.23489881f, 0.57123565f,  0.80051137f,  -1.28411946f,
  -1.36576732f, -1.30257508f, -1.30575106f
};

static const float av1_var_part_nn_bias_32_layer0[8] = {
  -1.6301435f, 0.61879037f, -1.68612662f, 1.66960165f,
  -0.0838243f, 0.32253287f, -0.65755282f, 0.96661531f
};

static const float av1_var_part_nn_weights_32_layer1[8] = {
  1.99257161f,  0.7331492f,  1.33539961f,  1.13501456f,
  -2.21154528f, 1.85858542f, -0.85565298f, -1.96410246f
};

static const float av1_var_part_nn_bias_32_layer1[1] = { -0.14880827f };

static const float av1_var_part_means_32[FEATURES] = {
  5.36360686f, 9.88421868f, 0.23543671f, 0.23621205f, 0.23409667f, 0.22855539f
};

static const float av1_var_part_vars_32[FEATURES] = {
  0.89077225f, 2.32312894f, 0.02167654f, 0.02392842f, 0.02466495f, 0.02047641f
};
*/
/* 6 features - full model ----
static const float av1_var_part_nn_weights_32_layer0[FEATURES * 8] = {
  1.2332769f, -1.9463526f, 2.74521738f, 2.74163489f, 2.75303344f, 2.845425f,
  0.57715147f, -1.81804495f, -0.73064292f, 0.71112037f, 0.64908075f,
-1.68229805f, -1.56201086f, 3.27755738f, 1.79604989f, 1.76783558f,
1.98824437f, 2.07598388f, 0.56998521f, -1.45529645f, 3.59854093f,
3.67199736f, 3.81625344f, 3.72041698f, -0.09884246f, 0.50324735f, -0.45622227f,
0.27361789f, -1.16606044f, 1.54897116f, -1.35899018f, 2.67418748f,
2.55753111f, 2.36084638f, 2.80962405f, 2.74925323f, 1.5669641f,
-0.4444419f, 1.64530465f, 1.37103728f, 1.77357033f, 1.5668947f, 0.73836879f,
-0.99320512f, -0.3616217f, 0.07951171f, -0.49994665f, 1.44515015f
};

static const float av1_var_part_nn_bias_32_layer0[8] = {
  -1.73359168f, -5.74840277f, -2.41337997f, -0.893716f,  0.30473623f,
-0.59958086f,  1.04164836f,  2.05473346f
};

static const float av1_var_part_nn_weights_32_layer1[8] = {
  1.75633289f, -1.00962976f, -0.93689681f, -1.6849214f, 0.42677414f, 0.9297025f,
-0.46860007f, -0.48764451f
};

static const float av1_var_part_nn_bias_32_layer1[1] = { 0.70485651f };

static const float av1_var_part_means_32[FEATURES] = {
  6.27425945f, 10.44423587f, 0.23728807f, 0.23624903f, 0.23620058f, 0.22890869f
};

static const float av1_var_part_vars_32[FEATURES] = {
  1.30119694f, 3.49359101f, 0.0222531f, 0.02299042f, 0.02474269f, 0.02303395f
};
*/
/* 4 features -full
static const float av1_var_part_nn_weights_32_layer0[FEATURES * 8] = {
  -2.01751919f, 1.12266665f, 0.04822318f, 0.04584756f,
  2.02224693f, -3.11683435f, -0.31638801f, -0.85380258f,
  -1.0547558f, -0.417637f, 0.49192043f, 0.13962656f,
  -0.86821727f, 2.53349589f, -0.3371422f, -0.50339667f,
  1.58077279f, -3.27086515f, 0.16343147f, -0.31348124f,
  0.14113951f, -0.62297956f, 1.32198277f, 0.53000123f,
  1.01122859f, -1.1893636f, 1.15732193f, 0.01561495f,
  0.34348803f, 2.04918574f, -0.58191134f, -0.48704956f
};

static const float av1_var_part_nn_bias_32_layer0[8] = {
  -1.78816796f, -0.59977052f,  0.96602767f,  0.89452586f,  1.34142104f,
-0.900506f, -0.05543005f, 0.5342975f
};

static const float av1_var_part_nn_weights_32_layer1[8] = {
  -0.51775707f, 1.24847285f, 0.45472642f, 0.77676123f, -1.18161113f,
-0.28649776f, -0.40197474f, -0.7127957f
};

static const float av1_var_part_nn_bias_32_layer1[1] = { 0.47861401 };

static const float av1_var_part_means_32[FEATURES] = {
  6.11553734f, 10.34953436f,  0.10542593f,  0.41790953f
};

static const float av1_var_part_vars_32[FEATURES] = {
  1.39164003f, 3.63359383f, 0.00493842f, 0.01805877f
};
*/

/* 4 features full - rdo --------
static const float av1_var_part_nn_weights_32_layer0[FEATURES * 8] = {
  0.21075759f, 1.02094863f, 0.19197288f, -0.09309328f,
  1.22900231f, -2.10517443f, 0.24150384f, -0.11170809f,
  -0.11973619f, 0.55370099f, -1.79906016f, -0.50037356f,
  0.66044361f, -1.79385794f, -0.48097574f, 0.69754214f,
  -0.63248407f, 1.40846731f, 0.37901071f, 1.25432744f,
  1.17890282f, -2.01296073f, -0.26822965f, -0.54413169f,
  -0.48939658f, 0.41487176f, -1.58283086f, -0.94150364f,
  0.20486748f, -1.47577182f, -0.73225218f, 0.75752268f
};

static const float av1_var_part_nn_bias_32_layer0[8] = {
  -0.87687603f, 1.17281391f, 1.06895638f, -1.90222494f, 0.44769699f,
-0.17592551f,  2.0036429f, -2.64856142f
};

static const float av1_var_part_nn_weights_32_layer1[8] = {
  0.18946134f, -0.75018394f, -0.61359973f, 0.6785535f, 0.49355979f, 0.72268511f,
0.79068719f, -1.02147634f
};

static const float av1_var_part_nn_bias_32_layer1[1] = { -1.30099246f };

static const float av1_var_part_means_32[FEATURES] = {
  6.19629874f, 10.28781195f, 0.1067754f, 0.41902138f
};

static const float av1_var_part_vars_32[FEATURES] = {
  1.42797355f, 3.6670224f, 0.00471344f, 0.01803315f
};
*/

static const float av1_var_part_nn_weights_32_layer0[FEATURES * 8] = {
  -0.29984873f, 1.13980896f,  0.52260552f,  2.08949251f, 0.24815907f,
  -1.93632679f, 0.16443032f,  -0.73081627f, 0.52303043f, -0.49295586f,
  2.15808381f,  0.45323922f,  -0.27362541f, 0.17863734f, -2.1878087f,
  0.13013841f,  -0.25871385f, 0.52390641f,  0.62869556f, 1.74243655f,
  -1.1662365f,  2.25428889f,  -0.10619761f, 0.2575906f,  -1.92629082f,
  2.28867462f,  0.01933312f,  0.56224773f,  1.21814666f, -1.44398102f,
  0.25578015f,  0.72862206f
};

static const float av1_var_part_nn_bias_32_layer0[8] = {
  1.63384504f, 0.04999198f, -1.21030761f, -3.03635579f,
  2.29331094f, -2.3567227f, 0.18191926f,  1.82475211f
};

static const float av1_var_part_nn_weights_32_layer1[8] = {
  0.96202868f,  -0.34866915f, -0.24882474f, -0.78452111f,
  -0.92075708f, -0.6521511f,  0.57053166f,  -0.61062313f
};

static const float av1_var_part_nn_bias_32_layer1[1] = { 0.41205638 };

static const float av1_var_part_means_32[FEATURES] = {
  5.92367594f, 10.16122563f, 0.11118023f, 0.4177506f
};

static const float av1_var_part_vars_32[FEATURES] = { 1.6181967f, 3.08201656f,
                                                      0.00496669f,
                                                      0.01738489f };

static const NN_CONFIG av1_var_part_nnconfig_32 = {
  FEATURES,  // num_inputs
  1,         // num_outputs
  1,         // num_hidden_layers
  {
      8,
  },  // num_hidden_nodes
  {
      av1_var_part_nn_weights_32_layer0,
      av1_var_part_nn_weights_32_layer1,
  },
  {
      av1_var_part_nn_bias_32_layer0,
      av1_var_part_nn_bias_32_layer1,
  },
};
/*
static const float av1_var_part_nn_weights_16_layer0[FEATURES * 8] = {
  0.45118305f,  -0.22068295f, 0.4604435f,   -0.1446326f,  -0.15765035f,
  0.42260198f,  -0.0945916f,  0.49544996f,  0.62781567f,  -0.41564372f,
  -0.39103292f, 0.44407624f,  0.48382613f,  -0.85424238f, -0.00961433f,
  0.25383582f,  0.14403897f,  0.00901859f,  -0.83201967f, -0.19323284f,
  0.59271213f,  0.69487457f,  0.6897112f,   0.62768521f,  0.9204492f,
  -1.42448347f, -0.16491054f, -0.10114424f, -0.1069687f,  -0.11289049f,
  0.26290832f,  -0.41850393f, 0.17239733f,  0.41770622f,  0.43725942f,
  0.19362467f,  -0.35955731f, -0.899446f,   0.49726389f,  0.66569571f,
  0.65893982f,  0.53199654f,  -0.1158694f,  -0.26472603f, 0.4155923f,
  0.15059544f,  0.09596755f,  0.26247133f
};

static const float av1_var_part_nn_bias_16_layer0[8] = {
  1.64486321f, -0.11851574f, 1.29322833f,  -0.61193136f,
  0.33027532f, 1.04197232f,  -0.80716674f, 0.88681233f
};

static const float av1_var_part_nn_weights_16_layer1[8] = {
  -1.02832118f, 0.72800106f, -0.42904783f, 1.44490586f,
  -1.03888227f, -0.9023916f, -1.51543102f, -0.43059521f
};

static const float av1_var_part_nn_bias_16_layer1[1] = { -0.85087946f };

static const float av1_var_part_means_16[FEATURES] = {
  5.32551326f, 8.218448f, 0.21954822f, 0.22808377f, 0.23019798f, 0.22320699f
};

static const float av1_var_part_vars_16[FEATURES] = { 0.86806032f, 2.39938956f,
                                                      0.01958579f, 0.02437927f,
                                                      0.02420755f, 0.0192003f };
*/
/* 4 features -full  balanced ---- */
/*
static const float av1_var_part_nn_weights_16_layer0[FEATURES * 8] = {
  0.82450546f, -1.6895598f, 0.21642135f, -0.14443887f,
  -0.18035923f, -0.55240796f, 0.34666463f, 0.16069096f,
  0.55461297f, -0.93533721f, 0.09891572f, 0.02133696f,
  -0.61931713f, -0.35128206f, -0.453986f, -0.21651371f,
  -0.55959882f, 1.38134328f, 0.5952496f, 0.25055191f,
  -1.04636164f, -0.08672306f, 0.08479656f, -0.04823812f,
  0.73720988f, -1.46627247f, 0.31105437f, 0.14069369f,
  0.01704834f, 0.58639624f, 0.51496669f, 0.03216761f
};

static const float av1_var_part_nn_bias_16_layer0[8] = {
  -0.28342379f, 0.40087014f, 0.58890165f, -0.57883941f, 0.25592934f,
0.11025883f, 0.43725984f, 1.21082125f
};

static const float av1_var_part_nn_weights_16_layer1[8] = {
  1.96841153f, 0.38476373f, -1.10453144f, 1.48283501f, 1.09304046f,
-0.86025845f, -1.39816501f, -0.90548448f
};

static const float av1_var_part_nn_bias_16_layer1[1] = { 1.40576256f };

static const float av1_var_part_means_16[FEATURES] = {
  6.21759176f, 9.15846021f, 0.09138057f, 0.40205081f
};

static const float av1_var_part_vars_16[FEATURES] = {
1.35855017e+00f, 5.00234036e+00f, 4.42717193e-03f, 1.97911485e-02f };
*/
/*
static const float av1_var_part_nn_weights_16_layer0[FEATURES * 8] = {
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f
};

static const float av1_var_part_nn_bias_16_layer0[8] = {
  0.0f, 0.0f, 0.0f, 0.0f,
  0.0f, 0.0f, 0.0f, 0.0f,
};

static const float av1_var_part_nn_weights_16_layer1[8] = {
  0.0f, 0.0f, 0.0f, 0.0f,   0.0f, 0.0f, 0.0f, 0.0f,
};

static const float av1_var_part_nn_bias_16_layer1[1] = { -1.0f };

static const float av1_var_part_means_16[FEATURES] = {
    0.0f, 0.0f, 0.0f, 0.0f,
};

static const float av1_var_part_vars_16[FEATURES] = {   1.0f, 1.0f, 1.0f, 1.0f,
 };
*/

static const float av1_var_part_nn_weights_16_layer0[FEATURES * 8] = {
  0.60218645f,  -0.95168278f, 0.51634686f,  1.21407959f, -0.359106f,
  0.55789991f,  -0.98935846f, 0.71498108f,  0.34282938f, -1.33121535f,
  0.09070362f,  -0.08645423f, -0.64268357f, 0.52120687f, -0.48585938f,
  -0.12821742f, -0.6786013f,  1.12591442f,  0.1664033f,  0.68337361f,
  -0.125435f,   0.11083573f,  0.77590665f,  1.51007686f, -0.24878526f,
  0.1745726f,   -1.57295658f, -0.78550166f, 1.34981018f, -2.24290581f,
  0.24371403f,  0.0307447f
};

static const float av1_var_part_nn_bias_16_layer0[8] = {
  1.28291967f, -1.60513075f, -1.56013152f, 1.31266831f,
  0.56477992f, 0.14924097f,  0.00690593f,  0.87213861f
};

static const float av1_var_part_nn_weights_16_layer1[8] = {
  -0.78753057f, -0.54292971f, 0.82204531f,  0.65523398f,
  0.45723168f,  0.63344418f,  -0.33035848f, -0.8755908f
};

static const float av1_var_part_nn_bias_16_layer1[1] = { 0.51096042f };

static const float av1_var_part_means_16[FEATURES] = { 5.86825445f, 9.27798299f,
                                                       0.0948979f,
                                                       0.41481675f };

static const float av1_var_part_vars_16[FEATURES] = { 1.57502631f, 3.72615647f,
                                                      0.00463681f,
                                                      0.01878185f };

static const NN_CONFIG av1_var_part_nnconfig_16 = {
  FEATURES,  // num_inputs
  1,         // num_outputs
  1,         // num_hidden_layers
  {
      8,
  },  // num_hidden_nodes
  {
      av1_var_part_nn_weights_16_layer0,
      av1_var_part_nn_weights_16_layer1,
  },
  {
      av1_var_part_nn_bias_16_layer0,
      av1_var_part_nn_bias_16_layer1,
  },
};

#undef FEATURES

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_AV1_ML_PARTITION_MODELS_H_
