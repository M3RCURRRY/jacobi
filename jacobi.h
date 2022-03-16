#ifndef NULL
#define NULL (void*)0
#endif

#define p_s "115792089237316195423570985008687907853269984665640564039457584007913111864739"
#define q_s "4824670384888174809315457708695329493868526062531865828358260708486391273721"
#define b_s "115792089237316195423570985008687907853269984665640564039457584007913111864737"

#define x_s "86773867670483746082709700497332781598034016226947998080571074369332058445893"
#define y_s "112365783073842736181650671369471258524588261898469739668735228897200995462257"
#define t_s "34969497036250894576602400300903320615522974325083826128152497107656115683806"
#define z_s "53269948314907700562817814682691121416875131851387647397949804668930162464909"

#ifndef JACOBI_H
#define JACOBI_H

#include <openssl/bn.h>

struct Point{
    BIGNUM * x;
    BIGNUM * y;
    BIGNUM * t;
    BIGNUM * z;
};

struct Curve{
    BIGNUM * p;
    BIGNUM * q;
    BIGNUM * b;
};

void init_point_wpar(struct Point *);

void init_point(struct Point *, char *, char *, char *, char *);

void init_curve(struct Curve *);

struct Point add_points(struct Point, struct Point);

struct Point double_point(struct Point p);

struct Point multiple_point(struct Point p, BIGNUM * k);

int if_contains(struct Point P);

void free_curve(struct Curve *);

void free_point(struct Point *);

#endif
