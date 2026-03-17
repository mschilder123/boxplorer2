#include "hmd.h"

#if !defined(_WIN32)

#include "utils.h"

HMD *CreateHMD() {
  DEBUG("No HMD support atm.");
  return nullptr;
}

#else // _WIN32

// Some bits from OVR tinyroom sample
#include <assert.h>
#include <stdio.h>

#include "shader_procs.h"
#include <Extras/OVR_Math.h>
#include <OVR_CAPI_GL.h>
#include <SDL_opengl.h>

#include "utils.h"

#pragma comment(lib, "libovr.lib")

using namespace OVR;

struct OculusTextureBuffer {
  ovrSession Session;
  ovrTextureSwapChain ColorTextureChain;
  ovrTextureSwapChain DepthTextureChain;
  GLuint fboId;
  Sizei texSize;

  OculusTextureBuffer(ovrSession session, Sizei size, int sampleCount)
      : Session(session), ColorTextureChain(nullptr),
        DepthTextureChain(nullptr), fboId(0), texSize(0, 0) {
    assert(sampleCount <=
           1); // The code doesn't currently handle MSAA textures.

    texSize = size;

    // This texture isn't necessarily going to be a rendertarget,
    // but it usually is.
    assert(session); // No HMD? A little odd.

    ovrTextureSwapChainDesc desc = {};
    desc.Type = ovrTexture_2D;
    desc.ArraySize = 1;
    desc.Width = size.w;
    desc.Height = size.h;
    desc.MipLevels = 1;
    desc.Format = OVR_FORMAT_R8G8B8A8_UNORM_SRGB;
    desc.SampleCount = sampleCount;
    desc.StaticImage = ovrFalse;

    {
      ovrResult result =
          ovr_CreateTextureSwapChainGL(Session, &desc, &ColorTextureChain);

      int length = 0;
      ovr_GetTextureSwapChainLength(session, ColorTextureChain, &length);

      if (OVR_SUCCESS(result)) {
        for (int i = 0; i < length; ++i) {
          GLuint chainTexId;
          ovr_GetTextureSwapChainBufferGL(Session, ColorTextureChain, i,
                                          &chainTexId);
          glBindTexture(GL_TEXTURE_2D, chainTexId);

          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        }
      }
    }

    desc.Format = OVR_FORMAT_D32_FLOAT;

    {
      ovrResult result =
          ovr_CreateTextureSwapChainGL(Session, &desc, &DepthTextureChain);

      int length = 0;
      ovr_GetTextureSwapChainLength(session, DepthTextureChain, &length);

      if (OVR_SUCCESS(result)) {
        for (int i = 0; i < length; ++i) {
          GLuint chainTexId;
          ovr_GetTextureSwapChainBufferGL(Session, DepthTextureChain, i,
                                          &chainTexId);
          glBindTexture(GL_TEXTURE_2D, chainTexId);

          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        }
      }
    }

    glGenFramebuffers(1, &fboId);
  }

  ~OculusTextureBuffer() {
    if (ColorTextureChain) {
      ovr_DestroyTextureSwapChain(Session, ColorTextureChain);
      ColorTextureChain = nullptr;
    }
    if (DepthTextureChain) {
      ovr_DestroyTextureSwapChain(Session, DepthTextureChain);
      DepthTextureChain = nullptr;
    }
    if (fboId) {
      glDeleteFramebuffers(1, &fboId);
      fboId = 0;
    }
  }

  Sizei GetSize() const { return texSize; }

  void SetAndClearRenderSurface() {
    GLuint curColorTexId;
    GLuint curDepthTexId;
    {
      int curIndex;
      ovr_GetTextureSwapChainCurrentIndex(Session, ColorTextureChain,
                                          &curIndex);
      ovr_GetTextureSwapChainBufferGL(Session, ColorTextureChain, curIndex,
                                      &curColorTexId);
    }
    {
      int curIndex;
      ovr_GetTextureSwapChainCurrentIndex(Session, DepthTextureChain,
                                          &curIndex);
      ovr_GetTextureSwapChainBufferGL(Session, DepthTextureChain, curIndex,
                                      &curDepthTexId);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, fboId);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           curColorTexId, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                           curDepthTexId, 0);

    // glViewport(0, 0, texSize.w, texSize.h);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // glEnable(GL_FRAMEBUFFER_SRGB);

    // glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fboId);

    // CHECK_ERROR;
    // CHECK_FRAMEBUFFER;
  }

  void UnsetRenderSurface() {
    glBindFramebuffer(GL_FRAMEBUFFER, fboId);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                           0, 0);
  }

  void Commit() {
    ovr_CommitTextureSwapChain(Session, ColorTextureChain);
    ovr_CommitTextureSwapChain(Session, DepthTextureChain);
  }
};

class Quest3 : public HMD {
private:
  ovrSession session;
  ovrHmdDesc hmdDesc;
  long long frameIndex;
  OculusTextureBuffer *eyeRenderTexture[2];
  double sensorSampleTime;
  ovrEyeRenderDesc eyeRenderDesc[2];
  ovrPosef EyeRenderPose[2];
  ovrTrackingState trackingState;
  ovrGraphicsLuid luid;
  ovrInputState inputState;

  // button press/report tracking.
  unsigned int buttons;
  unsigned int touches;

public:
  bool ok;

  Quest3()
      : frameIndex(0), eyeRenderTexture{nullptr, nullptr}, ok(false),
        buttons(0), touches(0) {
    ovrInitParams initParams = {ovrInit_RequestVersion | ovrInit_FocusAware,
                                OVR_MINOR_VERSION, NULL, 0, 0};
    ovrResult result = ovr_Initialize(&initParams);
    if (OVR_SUCCESS(result)) {
      result = ovr_Create(&session, &luid);
      if (OVR_SUCCESS(result)) {
        hmdDesc = ovr_GetHmdDesc(session);
        ovr_SetTrackingOriginType(session, ovrTrackingOrigin_EyeLevel);
        ok = true;
      } else {
        DEBUG("Failed to create session.");
      }
    } else {
      DEBUG("Failed to initialize libOVR.");
    }
  }

  bool AllocateTextures() {
    for (int eye = 0; eye < 2; ++eye) {
      ovrSizei idealTextureSize = ovr_GetFovTextureSize(
          session, ovrEyeType(eye), hmdDesc.DefaultEyeFov[eye], 1);
      eyeRenderTexture[eye] =
          new OculusTextureBuffer(session, idealTextureSize, 1);

      if (!eyeRenderTexture[eye]->ColorTextureChain ||
          !eyeRenderTexture[eye]->DepthTextureChain) {
        DEBUG("Failed to create textures.");
        return false;
      }
    }
    return true;
  }

  void FreeTextures() {
    for (int eye = 0; eye < 2; ++eye) {
      delete eyeRenderTexture[eye];
      eyeRenderTexture[eye] = nullptr;
    }
  }

  // Single src sidebyside framebuffer, left/right
  bool BlitFrom(GLuint srcFbo, int width, int height) {
    ovrTimewarpProjectionDesc posTimewarpProjectionDesc = {};

    // Blit Scene to Eye Buffers
    for (int eye = 0; eye < 2; ++eye) {
      // Switch to eye render target
      eyeRenderTexture[eye]->SetAndClearRenderSurface();

      // Get (approximate) projection matrix (shader used different one..)
      Matrix4f proj = ovrMatrix4f_Projection(hmdDesc.DefaultEyeFov[eye], 0.2f,
                                             1000.0f, ovrProjection_None);
      posTimewarpProjectionDesc =
          ovrTimewarpProjectionDesc_FromProjection(proj, ovrProjection_None);

      // Blit eye halves
      glBindFramebuffer(GL_READ_FRAMEBUFFER, srcFbo);
      Sizei dstSize = eyeRenderTexture[eye]->GetSize();
      int h = height;
      int w = width / 2;

      // DEBUG("blit(0,0,%d,%d,0,0,%d,%d)", w, h, dstSize.w, dstSize.h);

      if (eye == 0) {
        glBlitFramebuffer(
            0, 0, w, h, 0, 0, dstSize.w, dstSize.h,
            GL_COLOR_BUFFER_BIT /*| GL_DEPTH_BUFFER_BIT, GL_NEAREST*/,
            GL_LINEAR);
      } else {
        glBlitFramebuffer(
            w, 0, w + w, h, 0, 0, dstSize.w, dstSize.h,
            GL_COLOR_BUFFER_BIT /*| GL_DEPTH_BUFFER_BIT, GL_NEAREST*/,
            GL_LINEAR);
      }

      // CHECK_ERROR;
      // CHECK_FRAMEBUFFER;

      // Avoids an error when calling SetAndClearRenderSurface during next
      // iteration. Without this, during the next while loop iteration
      // SetAndClearRenderSurface would bind a framebuffer with an invalid
      // COLOR_ATTACHMENT0 because the texture ID associated with
      // COLOR_ATTACHMENT0 had been unlocked by calling wglDXUnlockObjectsNV.
      eyeRenderTexture[eye]->UnsetRenderSurface();

      // Commit changes to the textures so they get picked up frame
      eyeRenderTexture[eye]->Commit();
    }

    // Do distortion rendering, Present and flush/sync

    ovrLayerEyeFovDepth ld = {};
    ld.Header.Type = ovrLayerType_EyeFovDepth;
    ld.Header.Flags = ovrLayerFlag_TextureOriginAtBottomLeft; // Because OpenGL.
    ld.ProjectionDesc = posTimewarpProjectionDesc;
    ld.SensorSampleTime = sensorSampleTime;

    for (int eye = 0; eye < 2; ++eye) {
      ld.ColorTexture[eye] = eyeRenderTexture[eye]->ColorTextureChain;
      ld.DepthTexture[eye] = eyeRenderTexture[eye]->DepthTextureChain;
      ld.Viewport[eye] = Recti(eyeRenderTexture[eye]->GetSize());
      ld.Fov[eye] = hmdDesc.DefaultEyeFov[eye];
      ld.RenderPose[eye] = EyeRenderPose[eye];
    }

    ovrLayerHeader *layers = &ld.Header;
    ovrResult result =
        ovr_SubmitFrame(session, frameIndex, nullptr, &layers, 1);

    // CHECK_ERROR;

    // exit the rendering loop if submit returns an error, will retry on
    // ovrError_DisplayLost

    if (!OVR_SUCCESS(result)) {
      DEBUG("!OVR_SUCCESS");
      return false;
    }

    frameIndex++;

    return true;
  }

  bool GetHeadPose(float hmd_quat[4], float hmd_pos[3], float *ipd) {
    eyeRenderDesc[0] =
        ovr_GetRenderDesc(session, ovrEye_Left, hmdDesc.DefaultEyeFov[0]);
    eyeRenderDesc[1] =
        ovr_GetRenderDesc(session, ovrEye_Right, hmdDesc.DefaultEyeFov[1]);

    ovrPosef HmdToEyePose[2] = {eyeRenderDesc[0].HmdToEyePose,
                                eyeRenderDesc[1].HmdToEyePose};

    ovr_GetEyePoses(session, frameIndex, ovrTrue, HmdToEyePose, EyeRenderPose,
                    &sensorSampleTime);

    *ipd = (Vector3f(EyeRenderPose[0].Position) -
            Vector3f(EyeRenderPose[1].Position))
               .Length();

    ovrQuatf lq = EyeRenderPose[0].Orientation;

    hmd_quat[0] = lq.x;
    hmd_quat[1] = lq.y;
    hmd_quat[2] = lq.z;
    hmd_quat[3] = lq.w;

    hmd_pos[0] =
        (EyeRenderPose[0].Position.x + EyeRenderPose[1].Position.x) / 2;
    hmd_pos[1] =
        (EyeRenderPose[0].Position.y + EyeRenderPose[1].Position.y) / 2;
    hmd_pos[2] =
        (EyeRenderPose[0].Position.z + EyeRenderPose[1].Position.z) / 2;

    double frameTime = ovr_GetPredictedDisplayTime(session, 0);
    trackingState = ovr_GetTrackingState(session, frameTime, ovrTrue);

    return true;
  }

  bool _getHandPose(ovrPoseStatef pose, float hmd_quat[4], float hmd_pos[3]) {
    hmd_quat[0] = pose.ThePose.Orientation.x;
    hmd_quat[1] = pose.ThePose.Orientation.y;
    hmd_quat[2] = pose.ThePose.Orientation.z;
    hmd_quat[3] = pose.ThePose.Orientation.w;

    if (hmd_pos != nullptr) {
      hmd_pos[0] = pose.ThePose.Position.x;
      hmd_pos[1] = pose.ThePose.Position.y;
      hmd_pos[2] = pose.ThePose.Position.z;
    }

    return true;
  }

  bool GetLeftHandPose(float hmd_quat[4], float hmd_pos[3]) {
    return _getHandPose(trackingState.HandPoses[ovrHand_Left], hmd_quat,
                        hmd_pos);
  }

  bool GetRightHandPose(float hmd_quat[4], float hmd_pos[3]) {
    return _getHandPose(trackingState.HandPoses[ovrHand_Right], hmd_quat,
                        hmd_pos);
  }

  bool RefreshInputState(void) {
    return (OVR_SUCCESS(
        ovr_GetInputState(session, ovrControllerType_Touch, &inputState)));
  }

  float RightThumbstickX() { return inputState.Thumbstick[ovrHand_Right].x; }

  float RightThumbstickY() { return inputState.Thumbstick[ovrHand_Right].y; }

  float LeftThumbstickX() { return inputState.Thumbstick[ovrHand_Left].x; }

  float LeftThumbstickY() { return inputState.Thumbstick[ovrHand_Left].y; }

  float LeftIndexTrigger(void) { return inputState.IndexTrigger[ovrHand_Left]; }

  float LeftHandTrigger(void) { return inputState.HandTrigger[ovrHand_Left]; }

  float RightIndexTrigger(void) {
    return inputState.IndexTrigger[ovrHand_Right];
  }

  float RightHandTrigger(void) { return inputState.HandTrigger[ovrHand_Right]; }

  bool _button(unsigned int mask) {
    unsigned int cur = inputState.Buttons & mask;
    bool result = cur && !(buttons & mask);
    buttons &= ~mask;
    buttons |= cur;
    return result;
  }

  bool ButtonA(void) { return _button(ovrButton_A); }
  bool ButtonB(void) { return _button(ovrButton_B); }
  bool ButtonX(void) { return _button(ovrButton_X); }
  bool ButtonY(void) { return _button(ovrButton_Y); }
  bool ButtonEnter(void) { return _button(ovrButton_Enter); }
  bool RightThumbButton(void) { return _button(ovrButton_RThumb); }
  bool LeftThumbButton(void) { return _button(ovrButton_LThumb); }

#if 0
  Matrix4f lproj = ovrMatrix4f_Projection(hmdDesc.DefaultEyeFov[0], 0.2f, 1000.0f, ovrProjection_None);
  Matrix4f rproj = ovrMatrix4f_Projection(hmdDesc.DefaultEyeFov[1], 0.2f, 1000.0f, ovrProjection_None);

  char buf[BUFSIZ];
  lproj.ToString(buf, sizeof(buf));
  DEBUG("lproj %s", buf);
  rproj.ToString(buf, sizeof(buf));
  DEBUG("rproj %s", buf);
#endif

  ~Quest3() {
    if (session) {
      FreeTextures();
      ovr_Destroy(session);
      session = nullptr;
    }
    ovr_Shutdown();
  }

  void ResetPose() { ovr_RecenterTrackingOrigin(session); }
};

// Static factory
HMD *CreateHMD() {
  Quest3 *quest3 = new Quest3();
  if (!quest3->ok) {
    delete quest3;
    quest3 = nullptr;
  }
  return quest3;
}

#endif // _WIN32
