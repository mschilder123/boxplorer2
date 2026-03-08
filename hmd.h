#pragma once

class HMD {
public:
  virtual ~HMD() {}

  virtual bool AllocateTextures() = 0;
  virtual void FreeTextures() = 0;

  virtual bool BlitFrom(unsigned int srcFbo, int width, int height) = 0;

  virtual bool GetHeadPose(float hmd_quat[4], float hmd_pos[3]) = 0;
  virtual bool GetLeftHandPose(float hmd_quat[4], float hmd_pos[3]) = 0;
  virtual bool GetRightHandPose(float hmd_quat[4], float hmd_pos[3]) = 0;

  virtual void ResetPose() = 0;

  // controller acccess
  virtual bool RefreshInputState() = 0;

  virtual bool ButtonA() = 0;
  virtual bool ButtonB() = 0;
  virtual bool RightThumbButton() = 0;
  virtual bool ButtonX() = 0;
  virtual bool ButtonY() = 0;
  virtual bool LeftThumbButton() = 0;
  virtual bool ButtonEnter() = 0;

  virtual float RightIndexTrigger() = 0;
  virtual float RightHandTrigger() = 0;

  virtual float LeftIndexTrigger() = 0;
  virtual float LeftHandTrigger() = 0;

  virtual float RightThumbstickX() = 0;
  virtual float RightThumbstickY() = 0;

  virtual float LeftThumbstickX() = 0;
  virtual float LeftThumbstickY() = 0;
};

HMD *CreateHMD(void);
