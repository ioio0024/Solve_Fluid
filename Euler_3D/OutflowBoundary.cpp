#include "OutflowBoundary.h"

void OutflowBoundaryXLeft::operator ()(Int32 is, Int32 ie,
                                       Int32 js, Int32 je,
                                       Int32 ks, Int32 ke,
                                       FluidQuantity3DArray *sia,
                                       FluidQuantity3DArray *via
                                       ) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = sia->d(ie, j, k) - sia->d(ie - 1, j, k);
        sia->d(i, j, k) = sia->d(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = sia->u(ie, j, k) - sia->u(ie - 1, j, k);
        sia->u(i, j, k) = sia->u(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = sia->v(ie, j, k) - sia->v(ie - 1, j, k);
        sia->v(i, j, k) = sia->v(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = sia->w(ie, j, k) - sia->w(ie - 1, j, k);
        sia->w(i, j, k) = sia->w(ie - 1, j, k) + ds * (ie - 1 - i);
        ds = sia->p(ie, j, k) - sia->p(ie - 1, j, k);
        sia->p(i, j, k) = sia->p(ie - 1, j, k) + ds * (ie - 1 - i);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        dv = via->d(ie, j, k) - via->d(ie - 1, j, k);
        via->d(i, j, k) = via->d(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = via->u(ie, j, k) - via->u(ie - 1, j, k);
        via->u(i, j, k) = via->u(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = via->v(ie, j, k) - via->v(ie - 1, j, k);
        via->v(i, j, k) = via->v(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = via->w(ie, j, k) - via->w(ie - 1, j, k);
        via->w(i, j, k) = via->w(ie - 1, j, k) + dv * (ie - 1 - i);
        dv = via->p(ie, j, k) - via->p(ie - 1, j, k);
        via->p(i, j, k) = via->p(ie - 1, j, k) + dv * (ie - 1 - i);
      }
    }
  }
}

void OutflowBoundaryXRight::operator ()(Int32 is, Int32 ie,
                                        Int32 js, Int32 je,
                                        Int32 ks, Int32 ke,
                                        FluidQuantity3DArray *sia,
                                        FluidQuantity3DArray *via
                                        ) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = sia->d(is, j, k) - sia->d(is - 1, j, k);
        sia->d(i, j, k) = sia->d(is, j, k) + ds * (i - is);
        ds = sia->u(is, j, k) - sia->u(is - 1, j, k);
        sia->u(i, j, k) = sia->u(is, j, k) + ds * (i - is);
        ds = sia->v(is, j, k) - sia->v(is - 1, j, k);
        sia->v(i, j, k) = sia->v(is, j, k) + ds * (i - is);
        ds = sia->w(is, j, k) - sia->w(is - 1, j, k);
        sia->w(i, j, k) = sia->w(is, j, k) + ds * (i - is);
        ds = sia->p(is, j, k) - sia->p(is - 1, j, k);
        sia->p(i, j, k) = sia->p(is, j, k) + ds * (i - is);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        dv = via->d(is - 1, j, k) - via->d(is - 2, j, k);
        via->d(i, j, k) = via->d(is - 1, j, k) + dv * (i - is + 1);
        dv = via->u(is - 1, j, k) - via->u(is - 2, j, k);
        via->u(i, j, k) = via->u(is - 1, j, k) + dv * (i - is + 1);
        dv = via->v(is - 1, j, k) - via->v(is - 2, j, k);
        via->v(i, j, k) = via->v(is - 1, j, k) + dv * (i - is + 1);
        dv = via->w(is - 1, j, k) - via->w(is - 2, j, k);
        via->w(i, j, k) = via->w(is - 1, j, k) + dv * (i - is + 1);
        dv = via->p(is - 1, j, k) - via->p(is - 2, j, k);
        via->p(i, j, k) = via->p(is - 1, j, k) + dv * (i - is + 1);
      }
    }
  }
}

void OutflowBoundaryYLeft::operator ()(Int32 is, Int32 ie,
                                       Int32 js, Int32 je,
                                       Int32 ks, Int32 ke,
                                       FluidQuantity3DArray *sia,
                                       FluidQuantity3DArray *via
                                       ) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = sia->d(i, je, k) - sia->d(i, je - 1, k);
        sia->d(i, j, k) = sia->d(i, je - 1, k) + ds * (je - 1 - j);
        ds = sia->u(i, je, k) - sia->u(i, je - 1, k);
        sia->u(i, j, k) = sia->u(i, je - 1, k) + ds * (je - 1 - j);
        ds = sia->v(i, je, k) - sia->v(i, je - 1, k);
        sia->v(i, j, k) = sia->v(i, je - 1, k) + ds * (je - 1 - j);
        ds = sia->w(i, je, k) - sia->w(i, je - 1, k);
        sia->w(i, j, k) = sia->w(i, je - 1, k) + ds * (je - 1 - j);
        ds = sia->p(i, je, k) - sia->p(i, je - 1, k);
        sia->p(i, j, k) = sia->p(i, je - 1, k) + ds * (je - 1 - j);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        dv = via->d(i, je, k) - via->d(i, je - 1, k);
        via->d(i, j, k) = via->d(i, je - 1, k) + dv * (je - 1 - j);
        dv = via->u(i, je, k) - via->u(i, je - 1, k);
        via->u(i, j, k) = via->u(i, je - 1, k) + dv * (je - 1 - j);
        dv = via->v(i, je, k) - via->v(i, je - 1, k);
        via->v(i, j, k) = via->v(i, je - 1, k) + dv * (je - 1 - j);
        dv = via->w(i, je, k) - via->w(i, je - 1, k);
        via->w(i, j, k) = via->w(i, je - 1, k) + dv * (je - 1 - j);
        dv = via->p(i, je, k) - via->p(i, je - 1, k);
        via->p(i, j, k) = via->p(i, je - 1, k) + dv * (je - 1 - j);
      }
    }
  }
}

void OutflowBoundaryYRight::operator ()(Int32 is, Int32 ie,
                                        Int32 js, Int32 je,
                                        Int32 ks, Int32 ke,
                                        FluidQuantity3DArray *sia,
                                        FluidQuantity3DArray *via
                                        ) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = sia->d(i, js, k) - sia->d(i, js - 1, k);
        sia->d(i, j, k) = sia->d(i, js, k) + ds * (j - js);
        ds = sia->u(i, js, k) - sia->u(i, js - 1, k);
        sia->u(i, j, k) = sia->u(i, js, k) + ds * (j - js);
        ds = sia->v(i, js, k) - sia->v(i, js - 1, k);
        sia->v(i, j, k) = sia->v(i, js, k) + ds * (j - js);
        ds = sia->w(i, js, k) - sia->w(i, js - 1, k);
        sia->w(i, j, k) = sia->w(i, js, k) + ds * (j - js);
        ds = sia->p(i, js, k) - sia->p(i, js - 1, k);
        sia->p(i, j, k) = sia->p(i, js, k) + ds * (j - js);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        dv = via->d(i, js - 1, k) - via->d(i, js - 2, k);
        via->d(i, j, k) = via->d(i, js - 1, k) + dv * (j - js + 1);
        dv = via->u(i, js - 1, k) - via->u(i, js - 2, k);
        via->u(i, j, k) = via->u(i, js - 1, k) + dv * (j - js + 1);
        dv = via->v(i, js - 1, k) - via->v(i, js - 2, k);
        via->v(i, j, k) = via->v(i, js - 1, k) + dv * (j - js + 1);
        dv = via->w(i, js - 1, k) - via->w(i, js - 2, k);
        via->w(i, j, k) = via->w(i, js - 1, k) + dv * (j - js + 1);
        dv = via->p(i, js - 1, k) - via->p(i, js - 2, k);
        via->p(i, j, k) = via->p(i, js - 1, k) + dv * (j - js + 1);
      }
    }
  }
}

void OutflowBoundaryZLeft::operator ()(Int32 is, Int32 ie,
                                       Int32 js, Int32 je,
                                       Int32 ks, Int32 ke,
                                       FluidQuantity3DArray *sia,
                                       FluidQuantity3DArray *via
                                       ) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = sia->d(i, j, ke) - sia->d(i, j, ke - 1);
        sia->d(i, j, k) = sia->d(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = sia->u(i, j, ke) - sia->u(i, j, ke - 1);
        sia->u(i, j, k) = sia->u(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = sia->v(i, j, ke) - sia->v(i, j, ke - 1);
        sia->v(i, j, k) = sia->v(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = sia->w(i, j, ke) - sia->w(i, j, ke - 1);
        sia->w(i, j, k) = sia->w(i, j, ke - 1) + ds * (ke - 1 - k);
        ds = sia->p(i, j, ke) - sia->p(i, j, ke - 1);
        sia->p(i, j, k) = sia->p(i, j, ke - 1) + ds * (ke - 1 - k);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        dv = via->d(i, j, ke) - via->d(i, j, ke - 1);
        via->d(i, j, k) = via->d(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = via->u(i, j, ke) - via->u(i, j, ke - 1);
        via->u(i, j, k) = via->u(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = via->v(i, j, ke) - via->v(i, j, ke - 1);
        via->v(i, j, k) = via->v(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = via->w(i, j, ke) - via->w(i, j, ke - 1);
        via->w(i, j, k) = via->w(i, j, ke - 1) + dv * (ke - 1 - k);
        dv = via->p(i, j, ke) - via->p(i, j, ke - 1);
        via->p(i, j, k) = via->p(i, j, ke - 1) + dv * (ke - 1 - k);
      }
    }
  }
}

void OutflowBoundaryZRight::operator ()(Int32 is, Int32 ie,
                                        Int32 js, Int32 je,
                                        Int32 ks, Int32 ke,
                                        FluidQuantity3DArray *sia,
                                        FluidQuantity3DArray *via
                                        ) noexcept
{
  Real ds;
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        ds = sia->d(i, j, ks) - sia->d(i, j, ks - 1);
        sia->d(i, j, k) = sia->d(i, j, ks) + ds * (k - ks);
        ds = sia->u(i, j, ks) - sia->u(i, j, ks - 1);
        sia->u(i, j, k) = sia->u(i, j, ks) + ds * (k - ks);
        ds = sia->v(i, j, ks) - sia->v(i, j, ks - 1);
        sia->v(i, j, k) = sia->v(i, j, ks) + ds * (k - ks);
        ds = sia->w(i, j, ks) - sia->w(i, j, ks - 1);
        sia->w(i, j, k) = sia->w(i, j, ks) + ds * (k - ks);
        ds = sia->p(i, j, ks) - sia->p(i, j, ks - 1);
        sia->p(i, j, k) = sia->p(i, j, ks) + ds * (k - ks);
      }
    }
  }
  Real dv;
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        dv = via->d(i, j, ks - 1) - via->d(i, j, ks - 2);
        via->d(i, j, k) = via->d(i, j, ks - 1) + dv * (k - ks + 1);
        dv = via->u(i, j, ks - 1) - via->u(i, j, ks - 2);
        via->u(i, j, k) = via->u(i, j, ks - 1) + dv * (k - ks + 1);
        dv = via->v(i, j, ks - 1) - via->v(i, j, ks - 2);
        via->v(i, j, k) = via->v(i, j, ks - 1) + dv * (k - ks + 1);
        dv = via->w(i, j, ks - 1) - via->w(i, j, ks - 2);
        via->w(i, j, k) = via->w(i, j, ks - 1) + dv * (k - ks + 1);
        dv = via->p(i, j, ks - 1) - via->p(i, j, ks - 2);
        via->p(i, j, k) = via->p(i, j, ks - 1) + dv * (k - ks + 1);
      }
    }
  }
}
