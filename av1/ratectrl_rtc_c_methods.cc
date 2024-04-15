#include "av1/ratectrl_rtc.h"
namespace aom {
extern "C"
{

    void * createAV1Controller(const AV1RateControlRtcConfig &rc_cfg)
    {
        std::unique_ptr<AV1RateControlRTC> controller;
        controller = AV1RateControlRTC::Create(rc_cfg);
        return controller.release();
    }

    bool UpdateRateControl_AV1(void *controller, const AV1RateControlRtcConfig &rc_cfg)
    {
        return static_cast<AV1RateControlRTC *>(controller )->UpdateRateControl(rc_cfg);
    }

    int GetQP_AV1(void *controller)
    {
        return static_cast<AV1RateControlRTC *>(controller )->GetQP();
    }

    AV1LoopfilterLevel GetLoopfilterLevel_AV1(void *controller)
    {
        return static_cast<AV1RateControlRTC *>(controller )->GetLoopfilterLevel();
    }

    FrameDropDecision ComputeQP_AV1(void *controller, const AV1FrameParamsRTC &frame_params)
    {
        return static_cast<AV1RateControlRTC *>(controller )->ComputeQP(frame_params);
    }

    void PostEncodeUpdate_AV1(void *controller, uint64_t encoded_frame_size)
    {
        return static_cast<AV1RateControlRTC *>(controller )->PostEncodeUpdate(encoded_frame_size);
    }

    bool GetSegmentationData_AV1(void *controller, AV1SegmentationData *segmentation_data)
    {
        return static_cast<AV1RateControlRTC *>(controller )->GetSegmentationData(segmentation_data);
    }

    AV1CdefInfo GetCdefInfo_AV1(void *controller)
    {
        return static_cast<AV1RateControlRTC *>(controller )->GetCdefInfo();
    }

    AV1RateControlRtcConfig* create_av1_rate_control_config()
    {
        return new AV1RateControlRtcConfig();
    }

}//extern c
}  // namespace aom
