package hu.david.giczi.mvmxpert.wgs;

public class ToWGS implements WGS84, IUGG67 {

    private final double y_EOV;
    private final double x_EOV;

    private final double h_EOV;

    public ToWGS(double y_EOV, double x_EOV, double h_EOV) {
        this.y_EOV = y_EOV;
        this.x_EOV = x_EOV;
        this.h_EOV = h_EOV;
    }

    public static void main(String[] args) {

    }

}
