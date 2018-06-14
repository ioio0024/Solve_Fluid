#ifndef OUTFLOWBOUNDARY_H
#define OUTFLOWBOUNDARY_H

#include "global_definitions.h"
#include "Vector3D.h"

void OutflowBoundaryXLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept;

void OutflowBoundaryXRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept;

void OutflowBoundaryYLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept;

void OutflowBoundaryYRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept;

void OutflowBoundaryZLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept;

void OutflowBoundaryZRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept;
#endif // OUTFLOWBOUNDARY_H
