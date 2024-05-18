package hu.david.giczi.mvmxpert.wgs;

public class ToWGS implements WGS84, IUGG67 {

    private final double a = 6378160.0;
    private final double b = 6356774.516;
    private final double[][] MATRIX_A = new double[15][7];
    private final double[] MATRIX_l = new double[15];
    private final double y_EOV;
    private final double x_EOV;

    private final double h_EOV;

    public ToWGS(double y_EOV, double x_EOV, double h_EOV) {
        this.y_EOV = y_EOV;
        this.x_EOV = x_EOV;
        this.h_EOV = h_EOV;
        createMatrixA();
        createMatrix_l();
        runMatrixProcess();
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
        double[][] At_A = new double[7][7];
        int rowIndex = 0;
        int columnIndex = 0;
        double element;
        for (int i = 0; i < 7; i++) {
            element = 0.0;
            for (int j = 0; j < 15; j++) {
                element += ( At[i][j] * MATRIX_A[j][i] );
            }
            if( rowIndex == 7){
                rowIndex = 0;
                columnIndex++;
            }
            At_A[rowIndex++][columnIndex] = element;
        }
        double[][] inverse_At_A = new double[7][7];

    }
    private void createMatrixA(){
        MATRIX_A[0][0] = 1.0;
        MATRIX_A[0][1] = 0.0;
        MATRIX_A[0][2] = 0.0;
        MATRIX_A[0][3] = WGS84.point1_X;
        MATRIX_A[0][4] = 0.0;
        MATRIX_A[0][5] = - 1 * WGS84.point1_Z;
        MATRIX_A[0][6] = WGS84.point1_Y;

        MATRIX_A[1][0] = 0.0;
        MATRIX_A[1][1] = 1.0;
        MATRIX_A[1][2] = 0.0;
        MATRIX_A[1][3] = WGS84.point1_Y;
        MATRIX_A[1][4] = WGS84.point1_Z;
        MATRIX_A[1][5] = 0.0;
        MATRIX_A[1][6] = - 1 * WGS84.point1_X;

        MATRIX_A[2][0] = 0.0;
        MATRIX_A[2][1] = 0.0;
        MATRIX_A[2][2] = 1.0;
        MATRIX_A[2][3] = WGS84.point1_Z;
        MATRIX_A[2][4] = - 1 * WGS84.point1_Y;
        MATRIX_A[2][5] = WGS84.point1_X;
        MATRIX_A[2][6] = 0.0;

        MATRIX_A[3][0] = 1.0;
        MATRIX_A[3][1] = 0.0;
        MATRIX_A[3][2] = 0.0;
        MATRIX_A[3][3] = WGS84.point2_X;
        MATRIX_A[3][4] = 0.0;
        MATRIX_A[3][5] = - 1 * WGS84.point2_Z;
        MATRIX_A[3][6] = WGS84.point2_Y;

        MATRIX_A[4][0] = 0.0;
        MATRIX_A[4][1] = 1.0;
        MATRIX_A[4][2] = 0.0;
        MATRIX_A[4][3] = WGS84.point2_Y;
        MATRIX_A[4][4] = WGS84.point2_Z;
        MATRIX_A[4][5] = 0.0;
        MATRIX_A[4][6] = - 1 * WGS84.point2_X;

        MATRIX_A[5][0] = 0.0;
        MATRIX_A[5][1] = 0.0;
        MATRIX_A[5][2] = 1.0;
        MATRIX_A[5][3] = WGS84.point2_Z;
        MATRIX_A[5][4] = - 1 * WGS84.point2_Y;
        MATRIX_A[5][5] = WGS84.point2_X;
        MATRIX_A[5][6] = 0.0;

        MATRIX_A[6][0] = 1.0;
        MATRIX_A[6][1] = 0.0;
        MATRIX_A[6][2] = 0.0;
        MATRIX_A[6][3] = WGS84.point3_X;
        MATRIX_A[6][4] = 0.0;
        MATRIX_A[6][5] = - 1 * WGS84.point3_Z;
        MATRIX_A[6][6] = WGS84.point3_Y;

        MATRIX_A[7][0] = 0.0;
        MATRIX_A[7][1] = 1.0;
        MATRIX_A[7][2] = 0.0;
        MATRIX_A[7][3] = WGS84.point3_Y;
        MATRIX_A[7][4] = WGS84.point3_Z;
        MATRIX_A[7][5] = 0.0;
        MATRIX_A[7][6] = - 1 * WGS84.point3_X;

        MATRIX_A[8][0] = 0.0;
        MATRIX_A[8][1] = 0.0;
        MATRIX_A[8][2] = 1.0;
        MATRIX_A[8][3] = WGS84.point3_Z;
        MATRIX_A[8][4] = - 1 * WGS84.point3_Y;
        MATRIX_A[8][5] = WGS84.point3_X;
        MATRIX_A[8][6] = 0.0;

        MATRIX_A[9][0] = 1.0;
        MATRIX_A[9][1] = 0.0;
        MATRIX_A[9][2] = 0.0;
        MATRIX_A[9][3] = WGS84.point4_X;
        MATRIX_A[9][4] = 0.0;
        MATRIX_A[9][5] = - 1 * WGS84.point4_Z;
        MATRIX_A[9][6] = WGS84.point4_Y;

        MATRIX_A[10][0] = 0.0;
        MATRIX_A[10][1] = 1.0;
        MATRIX_A[10][2] = 0.0;
        MATRIX_A[10][3] = WGS84.point4_Y;
        MATRIX_A[10][4] = WGS84.point4_Z;
        MATRIX_A[10][5] = 0.0;
        MATRIX_A[10][6] = - 1 * WGS84.point4_X;

        MATRIX_A[11][0] = 0.0;
        MATRIX_A[11][1] = 0.0;
        MATRIX_A[11][2] = 1.0;
        MATRIX_A[11][3] = WGS84.point4_Z;
        MATRIX_A[11][4] = - 1 * WGS84.point4_Y;
        MATRIX_A[11][5] = WGS84.point4_X;
        MATRIX_A[11][6] = 0.0;

        MATRIX_A[12][0] = 1.0;
        MATRIX_A[12][1] = 0.0;
        MATRIX_A[12][2] = 0.0;
        MATRIX_A[12][3] = WGS84.point5_X;
        MATRIX_A[12][4] = 0.0;
        MATRIX_A[12][5] = - 1 * WGS84.point5_Z;
        MATRIX_A[12][6] = WGS84.point5_Y;

        MATRIX_A[13][0] = 0.0;
        MATRIX_A[13][1] = 1.0;
        MATRIX_A[13][2] = 0.0;
        MATRIX_A[13][3] = WGS84.point5_Y;
        MATRIX_A[13][4] = WGS84.point5_Z;
        MATRIX_A[13][5] = 0.0;
        MATRIX_A[13][6] = - 1 * WGS84.point5_X;

        MATRIX_A[14][0] = 0.0;
        MATRIX_A[14][1] = 0.0;
        MATRIX_A[14][2] = 1.0;
        MATRIX_A[14][3] = WGS84.point5_Z;
        MATRIX_A[14][4] = - 1 * WGS84.point5_Y;
        MATRIX_A[14][5] = WGS84.point5_X;
        MATRIX_A[14][6] = 0.0;
    }

    private void createMatrix_l(){
        MATRIX_l[0] = IUGG67.point1_X;
        MATRIX_l[1] = IUGG67.point1_Y;
        MATRIX_l[2] = IUGG67.point1_Z;
        MATRIX_l[3] = IUGG67.point2_X;
        MATRIX_l[4] = IUGG67.point2_Y;
        MATRIX_l[5] = IUGG67.point2_Z;
        MATRIX_l[6] = IUGG67.point3_X;
        MATRIX_l[7] = IUGG67.point3_Y;
        MATRIX_l[8] = IUGG67.point3_Z;
        MATRIX_l[9] = IUGG67.point4_X;
        MATRIX_l[10] = IUGG67.point4_Y;
        MATRIX_l[11] = IUGG67.point4_Z;
        MATRIX_l[12] = IUGG67.point5_X;
        MATRIX_l[13] = IUGG67.point5_Y;
        MATRIX_l[14] = IUGG67.point5_Z;
    }
    public static void main(String[] args) {

        ToWGS toWGS = new ToWGS(0.0,0.0,0.0);

    }

}
