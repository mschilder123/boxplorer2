#ifndef _F_OCULUSSDK_H_
#define _F_OCULUSSDK_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  unsigned int h_resolution;
  unsigned int v_resolution;
  float h_screen_size;
  float v_screen_size;
  float interpupillary_distance;
  float lens_separation_distance;
  float eye_to_screen_distance;
  float distortion_k[4];
  float chrom_abr[4];
} hmd_settings_t;

int InitOculusSDK();
void GetOculusView(float view[3]);
bool GetOculusQuat(float quat[4]);
void ReleaseOculusSDK();
void SetOculusPrediction(float time);

int GetOculusDeviceInfo(hmd_settings_t *hmd_settings);

void ResetOculusOrientation();

#ifdef __cplusplus
}
#endif

#endif
