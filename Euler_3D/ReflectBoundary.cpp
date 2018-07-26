#include "ReflectBoundary.h"

void ReflectBoundaryXLeft::operator ()(Int32 is, Int32 ie,
                                       Int32 js, Int32 je,
                                       Int32 ks, Int32 ke,
                                       FluidQuantity3DArray *sia,
                                       FluidQuantity3DArray *via
                                      ) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        sia->d(i, j, k) = sia->d(ie - 1, j, k);
        sia->u(i, j, k) = 0.0;
        sia->v(i, j, k) = sia->v(ie - 1, j, k);
        sia->w(i, j, k) = sia->w(ie - 1, j, k);
        sia->p(i, j, k) = sia->p(ie - 1, j, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        via->d(i, j, k) = via->d(ie - 1, j, k);
        via->u(i, j, k) = -via->u(ie - 1, j, k);
        via->v(i, j, k) = via->v(ie - 1, j, k);
        via->w(i, j, k) = via->w(ie - 1, j, k);
        via->p(i, j, k) = via->p(ie - 1, j, k);
      }
    }
  }
}

void ReflectBoundaryXRight::operator ()(Int32 is, Int32 ie,
                                        Int32 js, Int32 je,
                                        Int32 ks, Int32 ke,
                                        FluidQuantity3DArray *sia,
                                        FluidQuantity3DArray *via
                                       ) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        sia->d(i, j, k) = sia->d(is, j, k);
        sia->u(i, j, k) = 0.0;
        sia->v(i, j, k) = sia->v(is, j, k);
        sia->w(i, j, k) = sia->w(is, j, k);
        sia->p(i, j, k) = sia->p(is, j, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie - 1; ++i){
        via->d(i, j, k) = via->d(is - 1, j, k);
        via->u(i, j, k) = -via->u(is - 1, j, k);
        via->v(i, j, k) = via->v(is - 1, j, k);
        via->w(i, j, k) = via->w(is - 1, j, k);
        via->p(i, j, k) = via->p(is - 1, j, k);
      }
    }
  }
}

void ReflectBoundaryYLeft::operator ()(Int32 is, Int32 ie,
                                       Int32 js, Int32 je,
                                       Int32 ks, Int32 ke,
                                       FluidQuantity3DArray *sia,
                                       FluidQuantity3DArray *via
                                      ) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        sia->d(i, j, k) = sia->d(i, je - 1, k);
        sia->u(i, j, k) = sia->u(i, je - 1, k);
        sia->v(i, j, k) = 0.0;
        sia->w(i, j, k) = sia->w(i, je - 1, k);
        sia->p(i, j, k) = sia->p(i, je - 1, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        via->d(i, j, k) = via->d(i, je - 1, k);
        via->u(i, j, k) = via->u(i, je - 1, k);
        via->v(i, j, k) = -via->v(i, je - 1, k);
        via->w(i, j, k) = via->w(i, je - 1, k);
        via->p(i, j, k) = via->p(i, je - 1, k);
      }
    }
  }
}

void ReflectBoundaryYRight::operator ()(Int32 is, Int32 ie,
                                        Int32 js, Int32 je,
                                        Int32 ks, Int32 ke,
                                        FluidQuantity3DArray *sia,
                                        FluidQuantity3DArray *via
                                       ) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        sia->d(i, j, k) = sia->d(i, js, k);
        sia->u(i, j, k) = sia->u(i, js, k);
        sia->v(i, j, k) = 0.0;
        sia->w(i, j, k) = sia->w(i, js, k);
        sia->p(i, j, k) = sia->p(i, js, k);
      }
    }
  }
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je - 1; ++j){
      for(auto i = is; i < ie; ++i){
        via->d(i, j, k) = via->d(i, js - 1, k);
        via->u(i, j, k) = via->u(i, js - 1, k);
        via->v(i, j, k) = -via->v(i, js - 1, k);
        via->w(i, j, k) = via->w(i, js - 1, k);
        via->p(i, j, k) = via->p(i, js - 1, k);
      }
    }
  }
}

void ReflectBoundaryZLeft::operator ()(Int32 is, Int32 ie,
                                       Int32 js, Int32 je,
                                       Int32 ks, Int32 ke,
                                       FluidQuantity3DArray *sia,
                                       FluidQuantity3DArray *via
                                      ) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        sia->d(i, j, k) = sia->d(i, j, ke - 1);
        sia->u(i, j, k) = sia->u(i, j, ke - 1);
        sia->v(i, j, k) = sia->v(i, j, ke - 1);
        sia->w(i, j, k) = 0.0;
        sia->p(i, j, k) = sia->p(i, j, ke - 1);
      }
    }
  }
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        via->d(i, j, k) = via->d(i, j, ke - 1);
        via->u(i, j, k) = via->u(i, j, ke - 1);
        via->v(i, j, k) = via->v(i, j, ke - 1);
        via->w(i, j, k) = -via->w(i, j, ke - 1);
        via->p(i, j, k) = via->p(i, j, ke - 1);
      }
    }
  }
}

void ReflectBoundaryZRight::operator ()(Int32 is, Int32 ie,
                                        Int32 js, Int32 je,
                                        Int32 ks, Int32 ke,
                                        FluidQuantity3DArray *sia,
                                        FluidQuantity3DArray *via
                                       ) noexcept
{
  for(auto k = ks; k < ke; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        sia->d(i, j, k) = sia->d(i, j, ks);
        sia->u(i, j, k) = sia->u(i, j, ks);
        sia->v(i, j, k) = sia->v(i, j, ks);
        sia->w(i, j, k) = 0.0;
        sia->p(i, j, k) = sia->p(i, j, ks);
      }
    }
  }
  for(auto k = ks; k < ke - 1; ++k){
    for(auto j = js; j < je; ++j){
      for(auto i = is; i < ie; ++i){
        via->d(i, j, k) = via->d(i, j, ks - 1);
        via->u(i, j, k) = via->u(i, j, ks - 1);
        via->v(i, j, k) = via->v(i, j, ks - 1);
        via->w(i, j, k) = -via->w(i, j, ks - 1);
        via->p(i, j, k) = via->p(i, j, ks - 1);
      }
    }
  }
}
