#include <stdio.h>
#include "jacobi.h"

int main()
{
    struct Curve C = {NULL, NULL, NULL};
    struct Point P = {NULL, NULL, NULL, NULL};

    init_curve(&C);

    init_point_wpar(&P);

    //Tест 1
    BIGNUM * k = BN_new();
    BIGNUM * n = BN_new();
    BN_dec2bn(&n, "100000000000000000000000000");
    BN_rand_range(k,n);
    printf("1 тест. Является ли точка [k]P частью кривой: ");
    if (if_contains(multiple_point(P, k)) == 1) {
        printf("да\n");
    }
    else
        printf("нет\n");

    //Tест 2
    struct Point N = {NULL, NULL, NULL, NULL};
    init_point(&N, "0", "1", "1", "1");

    struct Point Test_2 = {NULL, NULL, NULL, NULL};
    Test_2 = multiple_point(P, C.q);

    printf("2 тест. Равны ли нейтральная точка и [q]P: ");
    if (BN_cmp(N.x, Test_2.x) && BN_cmp(N.y, Test_2.y) && BN_cmp(N.t, Test_2.t) && BN_cmp(N.z, Test_2.z))
        //if (N.x == Test_2.x && N.y == Test_2.y && N.t == Test_2.t && N.z == Test_2.z)
        printf("да\n");
    else
        printf("нет\n");

    //Tест 3
    printf("3.1 тест. Верно ли, что [q+1]P = P: ");
    BIGNUM * qq = BN_new();
    BN_dec2bn(&qq, "1");
    BN_add(qq, qq, C.q);
    struct Point Test_3 = {NULL, NULL, NULL, NULL};

    Test_3 = multiple_point(P, qq);

    if (BN_cmp(P.x,Test_3.x) && BN_cmp(P.y, Test_3.y) && BN_cmp(P.t, Test_3.t) && BN_cmp(P.z, Test_3.z))
        printf("да\n");
    else
        printf("нет\n");

    printf("3.2 тест. Верно ли, что [q-1]P = -P: ");
    BN_dec2bn(&qq, "1");
    BN_sub(qq, C.q, qq);

    Test_3 = multiple_point(P, qq);

    BN_dec2bn(&qq, "-1");
    BN_CTX *m_ctx = BN_CTX_new();
    BN_mul(qq, qq, P.x, m_ctx);
    if (BN_cmp(qq, Test_3.x) && BN_cmp(P.y, Test_3.y) && BN_cmp(P.t, Test_3.t) && BN_cmp(P.z, Test_3.z))
        printf("да\n");
    else
        printf("нет\n");

    //Tест 4
    BIGNUM * k1 = BN_new();
    BIGNUM * k2 = BN_new();
    BN_rand_range(k1, n);
    BN_rand_range(k2, n);
    struct Point Right = {NULL, NULL, NULL, NULL};
    struct Point Left = {NULL, NULL, NULL, NULL};

    Right = multiple_point(P, k2);
    Left = multiple_point(P, k1);
    Left = add_points(Left, Right);
    BN_add(k1, k1, k2);
    Right = multiple_point(P, k1);
    printf("4 тест. Верно ли, что [k1]P + [k2]P = [k1+k2]P: ");
    if (BN_cmp(Left.x, Right.x) && BN_cmp(Left.y, Right.y) && BN_cmp(Left.t, Right.t) && BN_cmp(Left.z, Right.z))
        printf("да\n");
    else
        printf("нет\n");

    BN_free(k);
    BN_free(n);
    BN_free(qq);
    BN_free(k1);
    BN_free(k2);

    BN_CTX_free(m_ctx);
}
