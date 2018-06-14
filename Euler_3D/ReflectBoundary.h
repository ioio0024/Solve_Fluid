#ifndef REFLECTBOUNDARY_H
#define REFLECTBOUNDARY_H

#include "global_definitions.h"
#include "Vector3D.h"

void ReflectBoundaryXLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept;

void ReflectBoundaryXRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept;

void ReflectBoundaryYLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept;

void ReflectBoundaryYRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept;

void ReflectBoundaryZLeft(Int32 is, Int32 ie,
                          Int32 js, Int32 je,
                          Int32 ks, Int32 ke,
                          Vector3DR *d_sia,
                          Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                          Vector3DR *p_sia,
                          Vector3DR *d_via,
                          Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                          Vector3DR *p_via) noexcept;

void ReflectBoundaryZRight(Int32 is, Int32 ie,
                           Int32 js, Int32 je,
                           Int32 ks, Int32 ke,
                           Vector3DR *d_sia,
                           Vector3DR *u_sia, Vector3DR *v_sia, Vector3DR *w_sia,
                           Vector3DR *p_sia,
                           Vector3DR *d_via,
                           Vector3DR *u_via, Vector3DR *v_via, Vector3DR *w_via,
                           Vector3DR *p_via) noexcept;
#endif // REFLECTBOUNDARY_H
