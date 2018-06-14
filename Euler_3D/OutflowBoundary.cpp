#include "OutflowBoundary.h"

void OutflowBoundaryXLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = (*d_sia)(ie, j, k) - (*d_sia)(ie - 1, j, k);
        (*d_sia)(i, j, k) = (*d_sia)(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = (*u_sia)(ie, j, k) - (*u_sia)(ie - 1, j, k);
        (*u_sia)(i, j, k) = (*u_sia)(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = (*v_sia)(ie, j, k) - (*v_sia)(ie - 1, j, k);
        (*v_sia)(i, j, k) = (*v_sia)(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = (*w_sia)(ie, j, k) - (*w_sia)(ie - 1, j, k);
        (*w_sia)(i, j, k) = (*w_sia)(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = (*p_sia)(ie, j, k) - (*p_sia)(ie - 1, j, k);
        (*p_sia)(i, j, k) = (*p_sia)(ie - 1, j, k) + ds * (ie - 1 - i);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        dv = (*d_via)(ie, j, k) - (*d_via)(ie - 1, j, k);
        (*d_via)(i, j, k) = (*d_via)(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = (*u_via)(ie, j, k) - (*u_via)(ie - 1, j, k);
        (*u_via)(i, j, k) = (*u_via)(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = (*v_via)(ie, j, k) - (*v_via)(ie - 1, j, k);
        (*v_via)(i, j, k) = (*v_via)(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = (*w_via)(ie, j, k) - (*w_via)(ie - 1, j, k);
        (*w_via)(i, j, k) = (*w_via)(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = (*p_via)(ie, j, k) - (*p_via)(ie - 1, j, k);
        (*p_via)(i, j, k) = (*p_via)(ie - 1, j, k) + dv * (ie - 1 - i);
      }
    }
  }
}

void OutflowBoundaryXRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = (*d_sia)(is, j, k) - (*d_sia)(is - 1, j, k);
        (*d_sia)(i, j, k) = (*d_sia)(is, j, k) + ds * (i - is);
        ds = (*u_sia)(is, j, k) - (*u_sia)(is - 1, j, k);
        (*u_sia)(i, j, k) = (*u_sia)(is, j, k) + ds * (i - is);
        ds = (*v_sia)(is, j, k) - (*v_sia)(is - 1, j, k);
        (*v_sia)(i, j, k) = (*v_sia)(is, j, k) + ds * (i - is);
        ds = (*w_sia)(is, j, k) - (*w_sia)(is - 1, j, k);
        (*w_sia)(i, j, k) = (*w_sia)(is, j, k) + ds * (i - is);
        ds = (*p_sia)(is, j, k) - (*p_sia)(is - 1, j, k);
        (*p_sia)(i, j, k) = (*p_sia)(is, j, k) + ds * (i - is);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        dv = (*d_via)(is - 1, j, k) - (*d_via)(is - 2, j, k);
        (*d_via)(i, j, k) = (*d_via)(is - 1, j, k) + dv * (i - is + 1);
        dv = (*u_via)(is - 1, j, k) - (*u_via)(is - 2, j, k);
        (*u_via)(i, j, k) = (*u_via)(is - 1, j, k) + dv * (i - is + 1);
        dv = (*v_via)(is - 1, j, k) - (*v_via)(is - 2, j, k);
        (*v_via)(i, j, k) = (*v_via)(is - 1, j, k) + dv * (i - is + 1);
        dv = (*w_via)(is - 1, j, k) - (*w_via)(is - 2, j, k);
        (*w_via)(i, j, k) = (*w_via)(is - 1, j, k) + dv * (i - is + 1);
        dv = (*p_via)(is - 1, j, k) - (*p_via)(is - 2, j, k);
        (*p_via)(i, j, k) = (*p_via)(is - 1, j, k) + dv * (i - is + 1);
      }
    }
  }
}

void OutflowBoundaryYLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = (*d_sia)(i, je, k) - (*d_sia)(i, je - 1, k);
        (*d_sia)(i, j, k) = (*d_sia)(i, je - 1, k) + ds * (je - 1 - j);
        ds = (*u_sia)(i, je, k) - (*u_sia)(i, je - 1, k);
        (*u_sia)(i, j, k) = (*u_sia)(i, je - 1, k) + ds * (je - 1 - j);
        ds = (*v_sia)(i, je, k) - (*v_sia)(i, je - 1, k);
        (*v_sia)(i, j, k) = (*v_sia)(i, je - 1, k) + ds * (je - 1 - j);
        ds = (*w_sia)(i, je, k) - (*w_sia)(i, je - 1, k);
        (*w_sia)(i, j, k) = (*w_sia)(i, je - 1, k) + ds * (je - 1 - j);
        ds = (*p_sia)(i, je, k) - (*p_sia)(i, je - 1, k);
        (*p_sia)(i, j, k) = (*p_sia)(i, je - 1, k) + ds * (je - 1 - j);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        dv = (*d_via)(i, je, k) - (*d_via)(i, je - 1, k);
        (*d_via)(i, j, k) = (*d_via)(i, je - 1, k) + dv * (je - 1 - j);
        dv = (*u_via)(i, je, k) - (*u_via)(i, je - 1, k);
        (*u_via)(i, j, k) = (*u_via)(i, je - 1, k) + dv * (je - 1 - j);
        dv = (*v_via)(i, je, k) - (*v_via)(i, je - 1, k);
        (*v_via)(i, j, k) = (*v_via)(i, je - 1, k) + dv * (je - 1 - j);
        dv = (*w_via)(i, je, k) - (*w_via)(i, je - 1, k);
        (*w_via)(i, j, k) = (*w_via)(i, je - 1, k) + dv * (je - 1 - j);
        dv = (*p_via)(i, je, k) - (*p_via)(i, je - 1, k);
        (*p_via)(i, j, k) = (*p_via)(i, je - 1, k) + dv * (je - 1 - j);
      }
    }
  }
}

void OutflowBoundaryYRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = (*d_sia)(i, js, k) - (*d_sia)(i, js - 1, k);
        (*d_sia)(i, j, k) = (*d_sia)(i, js, k) + ds * (j - js);
        ds = (*u_sia)(i, js, k) - (*u_sia)(i, js - 1, k);
        (*u_sia)(i, j, k) = (*u_sia)(i, js, k) + ds * (j - js);
        ds = (*v_sia)(i, js, k) - (*v_sia)(i, js - 1, k);
        (*v_sia)(i, j, k) = (*v_sia)(i, js, k) + ds * (j - js);
        ds = (*w_sia)(i, js, k) - (*w_sia)(i, js - 1, k);
        (*w_sia)(i, j, k) = (*w_sia)(i, js, k) + ds * (j - js);
        ds = (*p_sia)(i, js, k) - (*p_sia)(i, js - 1, k);
        (*p_sia)(i, j, k) = (*p_sia)(i, js, k) + ds * (j - js);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        dv = (*d_via)(i, js - 1, k) - (*d_via)(i, js - 2, k);
        (*d_via)(i, j, k) = (*d_via)(i, js - 1, k) + dv * (j - js + 1);
        dv = (*u_via)(i, js - 1, k) - (*u_via)(i, js - 2, k);
        (*u_via)(i, j, k) = (*u_via)(i, js - 1, k) + dv * (j - js + 1);
        dv = (*v_via)(i, js - 1, k) - (*v_via)(i, js - 2, k);
        (*v_via)(i, j, k) = (*v_via)(i, js - 1, k) + dv * (j - js + 1);
        dv = (*w_via)(i, js - 1, k) - (*w_via)(i, js - 2, k);
        (*w_via)(i, j, k) = (*w_via)(i, js - 1, k) + dv * (j - js + 1);
        dv = (*p_via)(i, js - 1, k) - (*p_via)(i, js - 2, k);
        (*p_via)(i, j, k) = (*p_via)(i, js - 1, k) + dv * (j - js + 1);
      }
    }
  }
}

void OutflowBoundaryZLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = (*d_sia)(i, j, ke) - (*d_sia)(i, j, ke - 1);
        (*d_sia)(i, j, k) = (*d_sia)(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = (*u_sia)(i, j, ke) - (*u_sia)(i, j, ke - 1);
        (*u_sia)(i, j, k) = (*u_sia)(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = (*v_sia)(i, j, ke) - (*v_sia)(i, j, ke - 1);
        (*v_sia)(i, j, k) = (*v_sia)(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = (*w_sia)(i, j, ke) - (*w_sia)(i, j, ke - 1);
        (*w_sia)(i, j, k) = (*w_sia)(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = (*p_sia)(i, j, ke) - (*p_sia)(i, j, ke - 1);
        (*p_sia)(i, j, k) = (*p_sia)(i, j, ke - 1) + ds * (ke - 1 - k);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        dv = (*d_via)(i, j, ke) - (*d_via)(i, j, ke - 1);
        (*d_via)(i, j, k) = (*d_via)(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = (*u_via)(i, j, ke) - (*u_via)(i, j, ke - 1);
        (*u_via)(i, j, k) = (*u_via)(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = (*v_via)(i, j, ke) - (*v_via)(i, j, ke - 1);
        (*v_via)(i, j, k) = (*v_via)(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = (*w_via)(i, j, ke) - (*w_via)(i, j, ke - 1);
        (*w_via)(i, j, k) = (*w_via)(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = (*p_via)(i, j, ke) - (*p_via)(i, j, ke - 1);
        (*p_via)(i, j, k) = (*p_via)(i, j, ke - 1) + dv * (ke - 1 - k);
      }
    }
  }
}

void OutflowBoundaryZRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = (*d_sia)(i, j, ks) - (*d_sia)(i, j, ks - 1);
        (*d_sia)(i, j, k) = (*d_sia)(i, j, ks) + ds * (k - ks);
        ds = (*u_sia)(i, j, ks) - (*u_sia)(i, j, ks - 1);
        (*u_sia)(i, j, k) = (*u_sia)(i, j, ks) + ds * (k - ks);
        ds = (*v_sia)(i, j, ks) - (*v_sia)(i, j, ks - 1);
        (*v_sia)(i, j, k) = (*v_sia)(i, j, ks) + ds * (k - ks);
        ds = (*w_sia)(i, j, ks) - (*w_sia)(i, j, ks - 1);
        (*w_sia)(i, j, k) = (*w_sia)(i, j, ks) + ds * (k - ks);
        ds = (*p_sia)(i, j, ks) - (*p_sia)(i, j, ks - 1);
        (*p_sia)(i, j, k) = (*p_sia)(i, j, ks) + ds * (k - ks);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        dv = (*d_via)(i, j, ks - 1) - (*d_via)(i, j, ks - 2);
        (*d_via)(i, j, k) = (*d_via)(i, j, ks - 1) + dv * (k - ks + 1);
        dv = (*u_via)(i, j, ks - 1) - (*u_via)(i, j, ks - 2);
        (*u_via)(i, j, k) = (*u_via)(i, j, ks - 1) + dv * (k - ks + 1);
        dv = (*v_via)(i, j, ks - 1) - (*v_via)(i, j, ks - 2);
        (*v_via)(i, j, k) = (*v_via)(i, j, ks - 1) + dv * (k - ks + 1);
        dv = (*w_via)(i, j, ks - 1) - (*w_via)(i, j, ks - 2);
        (*w_via)(i, j, k) = (*w_via)(i, j, ks - 1) + dv * (k - ks + 1);
        dv = (*p_via)(i, j, ks - 1) - (*p_via)(i, j, ks - 2);
        (*p_via)(i, j, k) = (*p_via)(i, j, ks - 1) + dv * (k - ks + 1);
      }
    }
  }
}
