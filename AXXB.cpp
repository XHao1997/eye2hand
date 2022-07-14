//
// Created by Hao on 25/06/2022.
//
# include "AXXB.h"

using namespace Eigen;
using namespace std;

Eigen::Matrix3d Eye2HandCli(Eigen::MatrixXd A, Eigen::MatrixXd B) {
    std::cout << "test0" << endl;
    int m = A.rows();
    int n = A.cols();
    cout << "m rows " << m << "n columns " << n << endl;
    Eigen::MatrixXd AA(m, n);
    Eigen::MatrixXd A1(3, 3);
    Eigen::MatrixXd B1(3, 3);
    int num = m / 4;
    for (int i = 0; i < num; i++) {
        A1 = A.block<3, 3>(4 * i, 0);
        B1 = B.block<3, 3>(4 * i, 0);

        Eigen::Quaterniond A2q = rotationMatrix2Quaterniond(A1);
        Eigen::Quaterniond B2q = rotationMatrix2Quaterniond(B1);
        AA.block<4, 4>(4 * i, 0) = CalAA(A2q, B2q);
    }
    cout << "AA: " << endl;
    cout << AA << endl;
    JacobiSVD<MatrixXd> svd(AA, ComputeFullU | ComputeFullV);
    MatrixXd V = svd.matrixV();
    Quaterniond v(V(0, 3), V(1, 3), V(2, 3), V(3, 3));
    cout << "Matrix V" << endl << V << endl;
    cout << "quaternion v" << endl << v << endl;
    Eigen::Matrix3d R = v.normalized().toRotationMatrix();
    cout << "Calculate the best R" << endl << R << endl;

    MatrixXd C(3 * num, 3);
    MatrixXd d(3 * num, 1);
    MatrixXd I(3, 3);
    I = MatrixXd::Identity(3, 3);
    for (int j = 0; j < num; j++) {
//        B1 = B.block<3,3>(4*j,0);
//        A1 = A.block<3,3>(4*j,0);
        C.block<3, 3>(3 * j, 0) = I - A.block<3, 3>(4 * j, 0);
        d.block<3, 1>(3 * j, 0) = A.block<3, 1>(4 * j, 3) - R * B.block<3, 1>(4 * j, 3);
    };
    MatrixXd t(3, 1);
    t = pinv_eigen_based(C, 0) * d;
    Matrix4d H;
    H.block<3, 3>(0, 0) = R;
    H.block<3, 1>(0, 3) = t;
    H.block<1, 4>(3, 0) << 0, 0, 0, 1;



//    std::cout<<"the size of C"<<endl<<C.rows()<<"  "<<C.cols()<<endl;
//    std::cout<<"the size of d"<<endl<<d.rows()<<"  "<<d.cols()<<endl;
    std::cout << "H is " << endl << H << endl;


    return R;
}

Eigen::Matrix4d CalAA(Eigen::Quaterniond a, Eigen::Quaterniond b) {
    Matrix4d RotMata;
    Matrix4d RotMatb;

    RotMata <<
            a.w(), -a.x(), -a.y(), -a.z(),
            a.x(), a.w() + skew(a)(0, 0), skew(a)(0, 1), skew(a)(0, 2),
            a.y(), skew(a)(1, 0), a.w() + skew(a)(1, 1), skew(a)(1, 2),
            a.z(), skew(a)(2, 0), skew(a)(2, 1), a.w() + skew(a)(2, 2);

    RotMatb <<
            b.w(), -b.x(), -b.y(), -b.z(),
            b.x(), b.w() - skew(b)(0, 0), -skew(b)(0, 1), -skew(b)(0, 2),
            b.y(), -skew(b)(1, 0), b.w() - skew(b)(1, 1), -skew(b)(1, 2),
            b.z(), -skew(b)(2, 0), -skew(b)(2, 1), b.w() - skew(b)(2, 2);
    Matrix4d c = RotMata - RotMatb;
    return c;
}

Eigen::Matrix3d skew(Eigen::Quaterniond V) {
    Eigen::Matrix3d skewV;
    skewV <<
          0, -V.z(), V.y(),
            V.z(), 0, -V.x(),
            -V.y(), V.x(), 0;
    return skewV;
}

Eigen::Quaterniond rotationMatrix2Quaterniond(Eigen::Matrix3d R) {
    Eigen::Quaterniond q = Eigen::Quaterniond(R);
    q.normalize();
    std::cout << "q is " << std::endl;
    std::cout << "x = " << q.x() << std::endl;
    std::cout << "y = " << q.y() << std::endl;
    std::cout << "z = " << q.z() << std::endl;
    std::cout << "w = " << q.w() << std::endl << std::endl;
    return q;
}


Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd &origin, const float er = 0) {
// 进行svd分解
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin, Eigen::ComputeThinU | Eigen::ComputeThinV);
// 构建SVD分解结果
    Eigen::MatrixXd U = svd_holder.matrixU();
    Eigen::MatrixXd V = svd_holder.matrixV();
    Eigen::MatrixXd D = svd_holder.singularValues();

// 构建S矩阵
    Eigen::MatrixXd S(V.cols(), U.cols());
    S.setZero();

    for (unsigned int i = 0; i < D.size(); ++i) {

        if (D(i, 0) > er) {
            S(i, i) = 1 / D(i, 0);
        } else {
            S(i, i) = 0;
        }
    }

// pinv_matrix = V * S * U^T
    return V * S * U.transpose();
}

int main() {
    Eigen::MatrixXd A(8, 4);
    Eigen::MatrixXd B(8, 4);
    Eigen::MatrixXd A2(4, 4);
    Eigen::MatrixXd B2(4, 4);
    Eigen::MatrixX4d A3(36, 4);
    Eigen::MatrixX4d B3(36, 4);
    A2 <<
       0.9916, 0.1018, -0.0800, 0.1046,
            -0.1018, 0.9948, 0.0046, -0.0298,
            0.0800, 0.0036, 0.9968, 0.0249,
            0, 0, 0, 1.000;
    B2 <<
       0.9926, 0.0698, -0.0988, -0.0629,
            -0.0652, 0.9966, 0.0495, -0.0770,
            0.1020, -0.0427, 0.9939, 0.1605,
            0, 0, 0, 1.000;
    A <<
      1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
    B <<
      1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
    cout << "--------------------------Test1---------------------" << endl;
    Eye2HandCli(A, B);
    cout << "--------------------------Test2---------------------" << endl;
    Eye2HandCli(A2, B2);
    cout << "--------------------------Test3---------------------" << endl;

    A3 << 0.991580492, 0.101840782, -0.07997864, 0.104560248,
            -0.101799113, 0.9947943, 0.004608923, -0.029840338,
            0.080031672, 0.003571636, 0.996785922, 0.024905401,
            0, 0, 0, 1,
            0.991803582, 0.097582025, -0.082482742, 0.09201491,
            -0.097390875, 0.995225969, 0.00634735, -0.044639178,
            0.082708354, 0.001737742, 0.99657228, 0.033178103,
            0, 0, 0, 1,
            0.992123279, 0.091537881, -0.0855115, 0.076112825,
            -0.09101402, 0.995799271, 0.010013012, -0.057967693,
            0.086068859, -0.002151397, 0.996286868, 0.038873647,
            0, 0, 0, 1,
            0.992537556, 0.084066533, -0.088329036, 0.057516068,
            -0.082967112, 0.996422999, 0.016051938, -0.069456555,
            0.089362514, -0.008603747, 0.995962006, 0.041756591,
            0, 0, 0, 1,
            0.993059607, 0.075585079, -0.090108337, 0.036936154,
            -0.073647188, 0.996979789, 0.024645319, -0.078932173,
            0.091699009, -0.017838045, 0.995626986, 0.041768556,
            0, 0, 0, 1,
            0.993715647, 0.066531362, -0.090015504, 0.015055714,
            -0.063540704, 0.997340711, 0.035694341, -0.086402994,
            0.09215092, -0.029750377, 0.995300519, 0.039035193,
            0, 0, 0, 1,
            0.994531963, 0.057328085, -0.087290692, -0.007538511,
            -0.05320689, 0.997388899, 0.048830439, -0.092018595,
            0.089862122, -0.043918966, 0.994985389, 0.033859199,
            0, 0, 0, 1,
            0.995514303, 0.048357312, -0.081319392, -0.030401127,
            -0.043248984, 0.997047575, 0.063448077, -0.096008566,
            0.084147481, -0.059646487, 0.994666526, 0.026699734,
            0, 0, 0, 1,
            0.996626628, 0.039950216, -0.071689225, -0.053245433,
            -0.034272234, 0.996304618, 0.078756096, -0.098611699,
            0.074570629, -0.076033472, 0.994312895, 0.018139972,
            0, 0, 0, 1;

    B3 <<
       0.992649553, 0.069827279, -0.098848453, -0.062860655,
            -0.065174127, 0.996643001, 0.049548589, -0.077019315,
            0.101976462, -0.042742023, 0.993868161, 0.160548521,
            0, 0, 0, 1,
            0.993052661, 0.06665936, -0.096968768, -0.083470284,
            -0.061716883, 0.996679717, 0.053109029, -0.062125057,
            0.100187018, -0.046755452, 0.993869453, 0.156594551,
            0, 0, 0, 1,
            0.993662328, 0.062819124, -0.093214462, -0.101734007,
            -0.057512036, 0.996625259, 0.05857012, -0.041197462,
            0.096579211, -0.052837968, 0.99392183, 0.143181675,
            0, 0, 0, 1,
            0.994461207, 0.058926538, -0.087032008, -0.11511078,
            -0.053299201, 0.996419789, 0.065626204, -0.0162654,
            0.09058754, -0.060623978, 0.994041564, 0.120051074,
            0, 0, 0, 1,
            0.995404846, 0.055533183, -0.078008059, -0.121596464,
            -0.04979272, 0.996036665, 0.073699711, 0.010170402,
            0.081791667, -0.069476816, 0.994224872, 0.087694206,
            0, 0, 0, 1,
            0.996411596, 0.053085881, -0.065922832, -0.119923213,
            -0.047612394, 0.995494986, 0.081992639, 0.035409181,
            0.0699785, -0.078559672, 0.994450294, 0.047317343,
            0, 0, 0, 1,
            0.997359339, 0.051911052, -0.050789671, -0.109665846,
            -0.047220242, 0.994861375, 0.089560558, 0.056842149,
            0.055177865, -0.086925759, 0.994685536, 0.000738204,
            0, 0, 0, 1,
            0.99809448, 0.052216633, -0.032876012, -0.091241964,
            -0.048872424, 0.994238256, 0.095403235, 0.07222668,
            0.037668225, -0.093614712, 0.994895668, -0.049771897,
            0, 0, 0, 1,
            0.998454492, 0.05410456, -0.012701324, -0.065809418,
            -0.052595264, 0.993739919, 0.098563241, 0.079906178,
            0.017954533, -0.097742881, 0.995049729, -0.101679436,
            0, 0, 0, 1;

    Eye2HandCli(A3, B3);
    // 设置矩阵行数、列数
    const int ROW = 3;
    const int COL = 4;

    // 生成大小 ROW * COL 的随机矩阵
    Eigen::MatrixXd A4;
    A4 = Eigen::MatrixXd::Random(ROW, COL);

    // 打印矩阵A
    cout << "origin is " << endl;
    cout << A4 << endl;

    // 打印矩阵A的伪逆矩阵
    cout << "pinv is " << endl;
    cout << pinv_eigen_based(A4) << endl;
    return 0;
}

