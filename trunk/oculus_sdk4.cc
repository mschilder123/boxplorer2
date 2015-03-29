// 2013 Dominic Szablewski - phoboslab.org
// sdk4 mods by marius

#include "oculus_sdk4.h"
#include "OVR.h"

#pragma comment(lib, "libovr.lib")
#pragma comment(lib, "winmm.lib")
#pragma comment(lib, "ole32.lib")

#ifndef M_PI
#define M_PI 3.14159265358979323846f // matches value in gcc v2 math.h
#endif

// Some libs that OVR needs..
#pragma comment(lib, "ws2_32.lib")
#pragma comment(lib, "gdi32.lib")

static ovrHmd hmd = NULL;

int InitOculusSDK() {
  ovr_Initialize();
  hmd = ovrHmd_Create(0);

  if (!hmd) {
    return 0;
  } else {
    ovrSizei resolution = hmd->Resolution;
    printf("HMD rez %dx%d\n", resolution.w, resolution.h);
  }

  ovrHmd_SetEnabledCaps(hmd,
      ovrHmdCap_LowPersistence | ovrHmdCap_DynamicPrediction);

  ovrHmd_ConfigureTracking(hmd,
      ovrTrackingCap_Orientation
      | ovrTrackingCap_MagYawCorrection
      | ovrTrackingCap_Position, 0);
  
  return 1;
}

void GetOculusView(float view[3]) {
#if 0
  // GetPredictedOrientation() works even if prediction is disabled
  OVR::Quatf q = fusion->GetPredictedOrientation();
  
  q.GetEulerAngles<OVR::Axis_Y, OVR::Axis_X, OVR::Axis_Z>(&view[1], &view[0], &view[2]);

  view[0] = (-view[0] * 180.0f) / M_PI;
  view[1] = (view[1] * 180.0f) / M_PI;
  view[2] = (-view[2] * 180.0f) / M_PI;
#endif
}

bool GetOculusQuat(float quat[4]) {
  bool result = false;
  ovrTrackingState ts = ovrHmd_GetTrackingState(hmd, 0.0/*ovr_GetTimeInSeconds()*/);
  if (ts.StatusFlags & (ovrStatus_OrientationTracked)) {
    ovrPoseStatef poseState = ts.HeadPose;
    ovrPosef posef = poseState.ThePose;
    ovrQuatf q = posef.Orientation;
  
    quat[0] = q.x;
    quat[1] = q.y;
    quat[2] = q.z;
    quat[3] = q.w;

    result = true;
  }
#if 1
  if (ts.StatusFlags & (ovrStatus_PositionTracked)) {
    ovrPoseStatef poseState = ts.HeadPose;
    ovrPosef posef = poseState.ThePose;
    ovrVector3f p = posef.Position;

    printf("pos(%f, %f, %f)\n", p.x, p.y, p.z);
  }
#endif
  return result;
}

void ReleaseOculusSDK()
{
  if (hmd) {
    ovrHmd_Destroy(hmd);
    hmd = NULL;
  }
  ovr_Shutdown();
}

void SetOculusPrediction(float time) {
}

int GetOculusDeviceInfo(hmd_settings_t *hmd_settings) {
  hmd_settings->h_resolution = hmd->Resolution.w;
  hmd_settings->v_resolution = hmd->Resolution.h;
#if 0
  hmd_settings->interpupillary_distance = hmdinfo.InterpupillaryDistance;
  hmd_settings->lens_separation_distance = hmdinfo.LensSeparationDistance;
  hmd_settings->eye_to_screen_distance = hmdinfo.EyeToScreenDistance;

  memcpy(hmd_settings->distortion_k,
                  hmdinfo.DistortionK, sizeof(float) * 4);
  memcpy(hmd_settings->chrom_abr,
                  hmdinfo.ChromaAbCorrection, sizeof(float) * 4);

#endif
  return 1;
}

void ResetOculusOrientation() {
  ovrHmd_RecenterPose(hmd);
}
