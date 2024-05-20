package hu.david.giczi.mvmxpert.wgs;

public class ToWGS implements WGS84, IUGG67 {

    private final double a = 6378137.0;
    private final double b = 6356752.314;
    public final double e = Math.sqrt((Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2));
    private final double e_ = Math.sqrt((Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(b, 2));
    private final double[][] MATRIX_A = new double[15][7];
    private final double[][] MATRIX_l = new double[15][1];
    private double[][] PARAM_FOR_WGS;
    public static double FI_WGS84;
    public static double LAMBDA_WGS84;
    public static double H_WGS84;


    public ToWGS(double X_IUGG67, double Y_IUGG67, double Z_IUGG67) {
        createMatrixA();
        createMatrix_l();
        runMatrixProcess();
        transformIUGG67CoordinatesForWGS(X_IUGG67, Y_IUGG67, Z_IUGG67);
    }

    private void transformIUGG67CoordinatesForWGS(double X_IUGG67, double Y_IUGG67, double Z_IUGG67){
      double deltaX = PARAM_FOR_WGS[0][0];
      double deltaY = PARAM_FOR_WGS[1][0];
      double deltaZ = PARAM_FOR_WGS[2][0];
      double k = PARAM_FOR_WGS[3][0];
      double eX = PARAM_FOR_WGS[4][0];
      double eY = PARAM_FOR_WGS[5][0];
      double eZ = PARAM_FOR_WGS[6][0];
      double[][] Rx =
                {{1.0, 0.0, 0.0},
                        {0.0, Math.cos(eX), Math.sin(eX)},
                        {0.0, - Math.sin(eX), Math.cos(eX)}};
      double[][] Ry =
                {{Math.cos(eY), 0.0, - Math.sin(eY)},
                        {0.0, 1.0, 0.0},
                        {Math.sin(eY), 0.0, Math.cos(eY)}};
      double[][] Rz =
                {{Math.cos(eZ), Math.sin(eZ), 0.0},
                        {- Math.sin(eZ), Math.cos(eZ), 0.0},
                        {0.0, 0.0, 1.0}};
        double kRxRy_00 = k * Rx[0][0] * Ry[0][0] + k * Rx[0][1] * Ry[1][0] + k * Rx[0][2] * Ry[2][0];
        double kRxRy_10 = k * Rx[1][0] * Ry[0][0] + k * Rx[1][1] * Ry[1][0] + k * Rx[1][2] * Ry[2][0];
        double kRxRy_20 = k * Rx[2][0] * Ry[0][0] + k * Rx[2][1] * Ry[1][0] + k * Rx[2][2] * Ry[2][0];

        double kRxRy_01 = k * Rx[0][0] * Ry[0][1] + k * Rx[0][1] * Ry[1][1] + k * Rx[0][2] * Ry[2][1];
        double kRxRy_11 = k * Rx[1][0] * Ry[0][1] + k * Rx[1][1] * Ry[1][1] + k * Rx[1][2] * Ry[2][1];
        double kRxRy_21 = k * Rx[2][0] * Ry[0][1] + k * Rx[2][1] * Ry[1][1] + k * Rx[2][2] * Ry[2][1];

        double kRxRy_02 = k * Rx[0][0] * Ry[0][2] + k * Rx[0][1] * Ry[1][2] + k * Rx[0][2] * Ry[2][2];
        double kRxRy_12 = k * Rx[1][0] * Ry[0][2] + k * Rx[1][1] * Ry[1][2] + k * Rx[1][2] * Ry[2][2];
        double kRxRy_22 = k * Rx[2][0] * Ry[0][2] + k * Rx[2][1] * Ry[1][2] + k * Rx[2][2] * Ry[2][2];

        double kRxRyRz_00 = kRxRy_00 * Rz[0][0] + kRxRy_01 * Rz[1][0] + kRxRy_02 * Rz[2][0];
        double kRxRyRz_10 = kRxRy_10 * Rz[0][0] + kRxRy_11 * Rz[1][0] + kRxRy_12 * Rz[2][0];
        double kRxRyRz_20 = kRxRy_20 * Rz[0][0] + kRxRy_21 * Rz[1][0] + kRxRy_22 * Rz[2][0];

        double kRxRyRz_01 = kRxRy_00 * Rz[0][1] + kRxRy_01 * Rz[1][1] + kRxRy_02 * Rz[2][1];
        double kRxRyRz_11 = kRxRy_10 * Rz[0][1] + kRxRy_11 * Rz[1][1] + kRxRy_12 * Rz[2][1];
        double kRxRyRz_21 = kRxRy_20 * Rz[0][1] + kRxRy_21 * Rz[1][1] + kRxRy_22 * Rz[2][1];

        double kRxRyRz_02 = kRxRy_00 * Rz[0][2] + kRxRy_01 * Rz[1][2] + kRxRy_02 * Rz[2][2];
        double kRxRyRz_12 = kRxRy_10 * Rz[0][2] + kRxRy_11 * Rz[1][2] + kRxRy_12 * Rz[2][2];
        double kRxRyRz_22 = kRxRy_20 * Rz[0][2] + kRxRy_21 * Rz[1][2] + kRxRy_22 * Rz[2][2];

        double X_WGS = deltaX + kRxRyRz_00 * X_IUGG67 + kRxRyRz_01 * Y_IUGG67 + kRxRyRz_02 * Z_IUGG67;
        double Y_WGS = deltaY + kRxRyRz_10 * X_IUGG67 + kRxRyRz_11 *  Y_IUGG67 + kRxRyRz_12 * Z_IUGG67;
        double Z_WGS = deltaZ + kRxRyRz_20 * X_IUGG67 + kRxRyRz_21 * Y_IUGG67 + kRxRyRz_22 * Z_IUGG67;

        double p = Math.sqrt(Math.pow(X_WGS, 2) + Math.pow(Y_WGS, 2));
        double theta = Math.atan(Z_WGS * a / (p * b));
        FI_WGS84 = Math.toDegrees(Math.atan( (Z_WGS + Math.pow(e_, 2) * b * Math.pow(Math.sin(theta), 3)) /
                (p - Math.pow(e, 2) * a * Math.pow(Math.cos(theta), 3)) ));
        LAMBDA_WGS84 = Math.toDegrees(Math.atan(Y_WGS / X_WGS));
        double N = a / Math.sqrt(1 - Math.pow(e, 2) * Math.pow(Math.sin(Math.toRadians(FI_WGS84)), 2));
        H_WGS84 = p / Math.cos(Math.toRadians(FI_WGS84)) - N;
    }

    private void runMatrixProcess() {
        //A transposed
        double[][] At = new double[7][15];
        for (int i = 0; i < MATRIX_A.length; i++) {
            for (int j = 0; j < MATRIX_A[i].length; j++) {
                At[j][i] = MATRIX_A[i][j];
            }
        }
        //At * A
        double[][] At_A = multiplyMatrix(At, MATRIX_A);
        double[][] inverse_At_A = InverseMatrix.invert(At_A);
        //At * l
        double[][] At_l = multiplyMatrix(At, MATRIX_l);
        //1 / ( At * A ) * (At * l)
        PARAM_FOR_WGS = multiplyMatrix(inverse_At_A, At_l);
    }

    private double[][] multiplyMatrix(double[][] matrix1, double[][] matrix2) {
        double[][] resultMatrix = new double[matrix1.length][matrix2[0].length];
        for (int i = 0; i < matrix1.length; i++) {
            for (int j = 0; j < matrix2[0].length; j++) {
                for (int k = 0; k < matrix2.length; k++) {
                    resultMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return resultMatrix;
    }

    private void createMatrixA(){
        MATRIX_A[0][0] = 1.0;
        MATRIX_A[0][1] = 0.0;
        MATRIX_A[0][2] = 0.0;
        MATRIX_A[0][3] = IUGG67.point1_X;
        MATRIX_A[0][4] = 0.0;
        MATRIX_A[0][5] = - 1 * IUGG67.point1_Z;
        MATRIX_A[0][6] = IUGG67.point1_Y;

        MATRIX_A[1][0] = 0.0;
        MATRIX_A[1][1] = 1.0;
        MATRIX_A[1][2] = 0.0;
        MATRIX_A[1][3] = IUGG67.point1_Y;
        MATRIX_A[1][4] = IUGG67.point1_Z;
        MATRIX_A[1][5] = 0.0;
        MATRIX_A[1][6] = - 1 * IUGG67.point1_X;

        MATRIX_A[2][0] = 0.0;
        MATRIX_A[2][1] = 0.0;
        MATRIX_A[2][2] = 1.0;
        MATRIX_A[2][3] = IUGG67.point1_Z;
        MATRIX_A[2][4] = - 1 * IUGG67.point1_Y;
        MATRIX_A[2][5] = IUGG67.point1_X;
        MATRIX_A[2][6] = 0.0;

        MATRIX_A[3][0] = 1.0;
        MATRIX_A[3][1] = 0.0;
        MATRIX_A[3][2] = 0.0;
        MATRIX_A[3][3] = IUGG67.point2_X;
        MATRIX_A[3][4] = 0.0;
        MATRIX_A[3][5] = - 1 * IUGG67.point2_Z;
        MATRIX_A[3][6] = IUGG67.point2_Y;

        MATRIX_A[4][0] = 0.0;
        MATRIX_A[4][1] = 1.0;
        MATRIX_A[4][2] = 0.0;
        MATRIX_A[4][3] = IUGG67.point2_Y;
        MATRIX_A[4][4] = IUGG67.point2_Z;
        MATRIX_A[4][5] = 0.0;
        MATRIX_A[4][6] = - 1 * IUGG67.point2_X;

        MATRIX_A[5][0] = 0.0;
        MATRIX_A[5][1] = 0.0;
        MATRIX_A[5][2] = 1.0;
        MATRIX_A[5][3] = IUGG67.point2_Z;
        MATRIX_A[5][4] = - 1 * IUGG67.point2_Y;
        MATRIX_A[5][5] = IUGG67.point2_X;
        MATRIX_A[5][6] = 0.0;

        MATRIX_A[6][0] = 1.0;
        MATRIX_A[6][1] = 0.0;
        MATRIX_A[6][2] = 0.0;
        MATRIX_A[6][3] = IUGG67.point3_X;
        MATRIX_A[6][4] = 0.0;
        MATRIX_A[6][5] = - 1 * IUGG67.point3_Z;
        MATRIX_A[6][6] = IUGG67.point3_Y;

        MATRIX_A[7][0] = 0.0;
        MATRIX_A[7][1] = 1.0;
        MATRIX_A[7][2] = 0.0;
        MATRIX_A[7][3] = IUGG67.point3_Y;
        MATRIX_A[7][4] = IUGG67.point3_Z;
        MATRIX_A[7][5] = 0.0;
        MATRIX_A[7][6] = - 1 * IUGG67.point3_X;

        MATRIX_A[8][0] = 0.0;
        MATRIX_A[8][1] = 0.0;
        MATRIX_A[8][2] = 1.0;
        MATRIX_A[8][3] = IUGG67.point3_Z;
        MATRIX_A[8][4] = - 1 * IUGG67.point3_Y;
        MATRIX_A[8][5] = IUGG67.point3_X;
        MATRIX_A[8][6] = 0.0;

        MATRIX_A[9][0] = 1.0;
        MATRIX_A[9][1] = 0.0;
        MATRIX_A[9][2] = 0.0;
        MATRIX_A[9][3] = IUGG67.point4_X;
        MATRIX_A[9][4] = 0.0;
        MATRIX_A[9][5] = - 1 * IUGG67.point4_Z;
        MATRIX_A[9][6] = IUGG67.point4_Y;

        MATRIX_A[10][0] = 0.0;
        MATRIX_A[10][1] = 1.0;
        MATRIX_A[10][2] = 0.0;
        MATRIX_A[10][3] = IUGG67.point4_Y;
        MATRIX_A[10][4] = IUGG67.point4_Z;
        MATRIX_A[10][5] = 0.0;
        MATRIX_A[10][6] = - 1 * IUGG67.point4_X;

        MATRIX_A[11][0] = 0.0;
        MATRIX_A[11][1] = 0.0;
        MATRIX_A[11][2] = 1.0;
        MATRIX_A[11][3] = IUGG67.point4_Z;
        MATRIX_A[11][4] = - 1 * IUGG67.point4_Y;
        MATRIX_A[11][5] = IUGG67.point4_X;
        MATRIX_A[11][6] = 0.0;

        MATRIX_A[12][0] = 1.0;
        MATRIX_A[12][1] = 0.0;
        MATRIX_A[12][2] = 0.0;
        MATRIX_A[12][3] = IUGG67.point5_X;
        MATRIX_A[12][4] = 0.0;
        MATRIX_A[12][5] = - 1 * IUGG67.point5_Z;
        MATRIX_A[12][6] = IUGG67.point5_Y;

        MATRIX_A[13][0] = 0.0;
        MATRIX_A[13][1] = 1.0;
        MATRIX_A[13][2] = 0.0;
        MATRIX_A[13][3] = IUGG67.point5_Y;
        MATRIX_A[13][4] = IUGG67.point5_Z;
        MATRIX_A[13][5] = 0.0;
        MATRIX_A[13][6] = - 1 * IUGG67.point5_X;

        MATRIX_A[14][0] = 0.0;
        MATRIX_A[14][1] = 0.0;
        MATRIX_A[14][2] = 1.0;
        MATRIX_A[14][3] = IUGG67.point5_Z;
        MATRIX_A[14][4] = - 1 * IUGG67.point5_Y;
        MATRIX_A[14][5] = IUGG67.point5_X;
        MATRIX_A[14][6] = 0.0;
    }

    private void createMatrix_l(){
        MATRIX_l[0][0] = WGS84.point1_X;
        MATRIX_l[1][0] = WGS84.point1_Y;
        MATRIX_l[2][0] = WGS84.point1_Z;
        MATRIX_l[3][0] = WGS84.point2_X;
        MATRIX_l[4][0] = WGS84.point2_Y;
        MATRIX_l[5][0] = WGS84.point2_Z;
        MATRIX_l[6][0] = WGS84.point3_X;
        MATRIX_l[7][0] = WGS84.point3_Y;
        MATRIX_l[8][0] = WGS84.point3_Z;
        MATRIX_l[9][0] = WGS84.point4_X;
        MATRIX_l[10][0] = WGS84.point4_Y;
        MATRIX_l[11][0] = WGS84.point4_Z;
        MATRIX_l[12][0] = WGS84.point5_X;
        MATRIX_l[13][0] = WGS84.point5_Y;
        MATRIX_l[14][0] = WGS84.point5_Z;
    }

}
