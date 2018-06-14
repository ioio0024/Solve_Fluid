#include "ReflectBoundary.h"

void ReflectBoundaryXLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_sia)(i, j, k) = (*d_sia)(ie - 1, j, k);
        (*u_sia)(i, j, k) = 0.0;
        (*v_sia)(i, j, k) = (*v_sia)(ie - 1, j, k);
        (*w_sia)(i, j, k) = (*w_sia)(ie - 1, j, k);
        (*p_sia)(i, j, k) = (*p_sia)(ie - 1, j, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        (*d_via)(i, j, k) = (*d_via)(ie - 1, j, k);
        (*u_via)(i, j, k) = -(*u_via)(ie - 1, j, k);
        (*v_via)(i, j, k) = (*v_via)(ie - 1, j, k);
        (*w_via)(i, j, k) = (*w_via)(ie - 1, j, k);
        (*p_via)(i, j, k) = (*p_via)(ie - 1, j, k);
      }
    }
  }
}

void ReflectBoundaryXRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_sia)(i, j, k) = (*d_sia)(is, j, k);
        (*u_sia)(i, j, k) = 0.0;
        (*v_sia)(i, j, k) = (*v_sia)(is, j, k);
        (*w_sia)(i, j, k) = (*w_sia)(is, j, k);
        (*p_sia)(i, j, k) = (*p_sia)(is, j, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        (*d_via)(i, j, k) = (*d_via)(is - 1, j, k);
        (*u_via)(i, j, k) = -(*u_via)(is - 1, j, k);
        (*v_via)(i, j, k) = (*v_via)(is - 1, j, k);
        (*w_via)(i, j, k) = (*w_via)(is - 1, j, k);
        (*p_via)(i, j, k) = (*p_via)(is - 1, j, k);
      }
    }
  }
}

void ReflectBoundaryYLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_sia)(i, j, k) = (*d_sia)(i, je - 1, k);
        (*u_sia)(i, j, k) = (*u_sia)(i, je - 1, k);
        (*v_sia)(i, j, k) = 0.0;
        (*w_sia)(i, j, k) = (*w_sia)(i, je - 1, k);
        (*p_sia)(i, j, k) = (*p_sia)(i, je - 1, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_via)(i, j, k) = (*d_via)(i, je - 1, k);
        (*u_via)(i, j, k) = (*u_via)(i, je - 1, k);
        (*v_via)(i, j, k) = -(*v_via)(i, je - 1, k);
        (*w_via)(i, j, k) = (*w_via)(i, je - 1, k);
        (*p_via)(i, j, k) = (*p_via)(i, je - 1, k);
      }
    }
  }
}

void ReflectBoundaryYRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_sia)(i, j, k) = (*d_sia)(i, js, k);
        (*u_sia)(i, j, k) = (*u_sia)(i, js, k);
        (*v_sia)(i, j, k) = 0.0;
        (*w_sia)(i, j, k) = (*w_sia)(i, js, k);
        (*p_sia)(i, j, k) = (*p_sia)(i, js, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_via)(i, j, k) = (*d_via)(i, js - 1, k);
        (*u_via)(i, j, k) = (*u_via)(i, js - 1, k);
        (*v_via)(i, j, k) = -(*v_via)(i, js - 1, k);
        (*w_via)(i, j, k) = (*w_via)(i, js - 1, k);
        (*p_via)(i, j, k) = (*p_via)(i, js - 1, k);
      }
    }
  }
}

void ReflectBoundaryZLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_sia)(i, j, k) = (*d_sia)(i, j, ke - 1);
        (*u_sia)(i, j, k) = (*u_sia)(i, j, ke - 1);
        (*v_sia)(i, j, k) = (*v_sia)(i, j, ke - 1);
        (*w_sia)(i, j, k) = 0.0;
        (*p_sia)(i, j, k) = (*p_sia)(i, j, ke - 1);
      }
    }
  }
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_via)(i, j, k) = (*d_via)(i, j, ke - 1);
        (*u_via)(i, j, k) = (*u_via)(i, j, ke - 1);
        (*v_via)(i, j, k) = (*v_via)(i, j, ke - 1);
        (*w_via)(i, j, k) = -(*w_via)(i, j, ke - 1);
        (*p_via)(i, j, k) = (*p_via)(i, j, ke - 1);
      }
    }
  }
}

void ReflectBoundaryZRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_sia)(i, j, k) = (*d_sia)(i, j, ks);
        (*u_sia)(i, j, k) = (*u_sia)(i, j, ks);
        (*v_sia)(i, j, k) = (*v_sia)(i, j, ks);
        (*w_sia)(i, j, k) = 0.0;
        (*p_sia)(i, j, k) = (*p_sia)(i, j, ks);
      }
    }
  }
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        (*d_via)(i, j, k) = (*d_via)(i, j, ks - 1);
        (*u_via)(i, j, k) = (*u_via)(i, j, ks - 1);
        (*v_via)(i, j, k) = (*v_via)(i, j, ks - 1);
        (*w_via)(i, j, k) = -(*w_via)(i, j, ks - 1);
        (*p_via)(i, j, k) = (*p_via)(i, j, ks - 1);
      }
    }
  }
}
