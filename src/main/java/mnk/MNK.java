package mnk;

import com.sun.tools.javac.util.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.knowm.xchart.Chart;
import org.knowm.xchart.SwingWrapper;

import java.lang.Double;
import java.util.*;
import java.util.List;

import org.knowm.xchart.StyleManager.LegendPosition;
import org.knowm.xchart.StyleManager.ChartType;

/**
 * Created by nikita on 18/10/15.
 */

public class MNK {

    static final int PERCENT = 20;

    //COSX
//    static final double TOP_Y = 1;
//    static final double BOTTOM_Y = -1;
//    static final double LEFT_X = 3.14;
//    static final double RIGHT_X = 3.14 * 4;
//    static final double RANGE_X = 0.25;
//    static final double MAX = 1.0;
//    static final int RANDOM_COUNT = 100;

    //    sin
//    static final double TOP_Y = 3.14;
//    static final double BOTTOM_Y = -3.14;
//    static final double LEFT_X = 0;
//    static final double RIGHT_X = 2.5;
//    static final double RANGE_X = 0.1;
//    static final double MAX = 3.14;
//    static final int RANDOM_COUNT = 90;

//    5x^3 + x^2 + 5
    static final double TOP_Y = 100;
    static final double BOTTOM_Y = -100;
    static final double LEFT_X = -7;
    static final double RIGHT_X = 7;
    static final double RANGE_X = 0.20;
    static final double MAX = 2000;
    static final int RANDOM_COUNT = 200;

    private final UniformRealDistribution distribution = new UniformRealDistribution(LEFT_X, RIGHT_X);
    private final NormalDistribution normalDistribution = new NormalDistribution();

    private List<List<Double>> funcData = solvData();
    private List<List<Double>> funcDataError = solvDataError();
    private int currentFunc = 3;


    private double[][] matrix;
    private int[] mask;

    double[][] makeSystem(double[][] xyTable, int basis) {

        double[][] matrix = new double[basis][basis+1];
        double lambda = 0.00001;

        for(int i = 0; i < basis; i++) {
            for(int j = 0; j < basis; j++) {
                double sumA = 0;
                double sumB = 0;
                for (int k = 0; k < xyTable[0].length; k++) {
                    sumA += Math.pow(xyTable[0][k], i) * Math.pow(xyTable[0][k], j);
                    sumB += xyTable[1][k] * Math.pow(xyTable[0][k], i);
                }
                matrix[i][j] = sumA;
                matrix[i][basis] = sumB;

                if (i == j) {
                    matrix[i][j] = sumA - lambda;
                }
            }
        }
        this.matrix = matrix;
        return matrix;
    }

    private double[] Gauss(int rowCount, int colCount){
        int i;
        mask = new int[colCount - 1];
        for (i = 0; i < colCount - 1; i++)
            mask[i] = i;
        if (GaussDirectPass(rowCount, colCount)) {
            double[] answer = GaussReversePass(colCount, rowCount);
            return answer;
        }
        else return null;
    }

    private boolean GaussDirectPass(int rowCount, int colCount){
        int i, j, k, maxId, tmpInt;
        double maxVal, tmpDouble;
        for (i = 0; i < rowCount; i++)
        {
            maxId = i;
            maxVal = matrix[i][i];
            for (j = i + 1; j < colCount - 1; j++)
                if(Math.abs(maxVal)<Math.abs(matrix[i][j]))
                {
                    maxVal = matrix[i][j];
                    maxId = j;
                }
            if (maxVal == 0) return false;
            if (i != maxId)
            {
                for (j = 0; j < rowCount; j++)
                {
                    tmpDouble = matrix[j][i];
                    matrix[j][i] = matrix[j][maxId];
                    matrix[j][maxId] = tmpDouble;
                }
                tmpInt = mask[i];
                mask[i] = mask[maxId];
                mask[maxId] = tmpInt;
            }
            for (j = 0; j < colCount; j++) matrix[i][j] /= maxVal;
            for (j = i + 1; j < rowCount; j++)
            {
                double tempMn = matrix[j][i];
                for (k = 0; k < colCount; k++)
                    matrix[j][k] -= matrix[i][k] * tempMn;
            }
        }
        return true;
    }

    private double[] GaussReversePass(int colCount, int rowCount){
        int i,j,k;
        for (i = rowCount - 1; i >= 0; i--)
            for (j = i - 1; j >= 0; j--)
            {
                double tempMn = matrix[j][i];
                for (k = 0; k < colCount; k++)
                    matrix[j][k] -= matrix[i][k] * tempMn;
            }
        double[] answer = new double[rowCount];
        for (i = 0; i < rowCount; i++) answer[mask[i]] = matrix[i][colCount - 1];
        return answer;
    }

    double[] LeastSquares(double[][] xyTable, int basis) {

        double[][] matrix = makeSystem(xyTable, basis);
//        double[] result = Kramer(matrix);
        double[] result = Gauss(basis, basis + 1);

        String str = "";
        for (int i = 0; i < basis-1; i++) {
            str += result[i] + "x^" + i + " + ";
        }
        str += result[basis-1] + "x^" + (basis-1);
//        System.out.println(str);
        return result;
    }

    double fy(double[]c, double x) {
        double sum = 0;
        for (int i = 0; i < c.length; i++) {
            sum += c[i] * Math.pow(x,i);
        }
        return sum;
    }

    double origFy(double x) {

        switch (currentFunc) {
            case 1: {
                return Math.cos(x);
            }
            case 2: {
                return  (x * Math.sin(2 * Math.PI * x));
            }
            case 3: {
                return (5*Math.pow(x,3) + Math.pow(x,2) + 5);
            }
        }

        return 0;
    }

    List<List<Double>> cosData() {

        List<Double> xData = new ArrayList<Double>();
        List<Double> yData = new ArrayList<Double>();

        for (double i = LEFT_X; i <= RIGHT_X; i+= RANGE_X) {
            xData.add(i);
            yData.add(Math.cos(i));
        }

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);

        return data;
    }

    List<List<Double>> cosDataError() {

        List<Double> xData = new ArrayList<Double>();
        List<Double> yData = new ArrayList<Double>();


        for (int i = 0; i < RANDOM_COUNT; i++) {
            double x = distribution.sample();
            double y = Math.cos(x);
            xData.add(x);
            yData.add(y);
        }

        for (int i = 0; i < yData.size(); i++) {
            if (new Random().nextInt(100) < PERCENT) {
                yData.set(i, normalDistribution.sample());
            }
        }

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);

        return data;
    }


    List<List<Double>> sinData() {

        List<Double> xData = new ArrayList<Double>();
        List<Double> yData = new ArrayList<Double>();

        for (double i = LEFT_X; i <= RIGHT_X; i+= RANGE_X) {
            xData.add(i);
            yData.add(i * Math.sin(2 * Math.PI * i));
        }

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);

        return data;
    }

    List<List<Double>> sinDataError() {

        List<Double> xData = new ArrayList<Double>();
        List<Double> yData = new ArrayList<Double>();

        for (int i = 0; i < RANDOM_COUNT; i++) {
            double x = distribution.sample();
            double y = x * Math.sin(2 * Math.PI * x);
            xData.add(x);
            yData.add(y);
        }

        for (int i = 0; i < yData.size(); i++) {
            if (new Random().nextInt(100) < PERCENT) {
                yData.set(i, normalDistribution.sample());
            }
        }

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);

        return data;
    }


    List<List<Double>> solvData() {

        List<Double> xData = new ArrayList<Double>();
        List<Double> yData = new ArrayList<Double>();

        for (double i = LEFT_X; i <= RIGHT_X; i+= RANGE_X) {
            xData.add(i);
            yData.add(5*Math.pow(i,3) + Math.pow(i,2) + 5);
        }

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);

        return data;
    }

    List<List<Double>> solvDataError() {

        List<Double> xData = new ArrayList<Double>();
        List<Double> yData = new ArrayList<Double>();

        for (int i = 0; i < RANDOM_COUNT; i++) {
            double x = distribution.sample();
            double y = 5 * Math.pow(x, 3) + Math.pow(x, 2) + 5;
            xData.add(x);
            yData.add(y);
        }

        for (int i = 0; i < yData.size(); i++) {
            if (new Random().nextInt(100) < PERCENT) {
                yData.set(i, normalDistribution.sample()*50);
            }
        }

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);

        return data;
    }

    public void crossValidation() {

        List<List<Double>> dataError = funcDataError;

        final int count = 8;
        int power = 3;

        double minError = Float.MAX_VALUE;
        double bestPower = 0;

        for (int i = 0; i < count; i++) {

            List<List<Double>> blocks = blocksWithoutIndex(count, i, dataError);

            double[] xa = xForList(blocks);
            double[] ya = yForList(blocks);
            double[][] xyTable = {xa, ya};
            double[] result = LeastSquares(xyTable, power);

            List<List<Double>> trainBlock = blockWithIndex(count, i, dataError);

            double[] trainX = xForList(trainBlock);

            double error = 0;
            for (int j = 0; j < trainX.length; j++) {
                double x = trainX[j];
                double y1 = fy(result, x);
                double y2 = origFy(x);
                error += (y1-y2)*(y1-y2);
            }
            error = error * 0.5f;

            if (minError > error) {
                minError = error;
                bestPower = power;
            }

            power++;
        }

        System.out.println("BEST POWER = " + bestPower);

    }

    public double[] yForList(List<List<Double>> list) {

        List<Double> yBlocks = list.get(1);
        Double[] yDataArray = yBlocks.toArray(new Double[yBlocks.size()]);
        double[] yDatadouble = new double[yDataArray.length];
        for (int i = 0; i < yDataArray.length; i++) {
            yDatadouble[i] = yDataArray[i];
        }
        return yDatadouble;
    }

    public double[] xForList(List<List<Double>> list) {

        List<Double> xBlocks = list.get(0);
        Double[] xDataArray = xBlocks.toArray(new Double[xBlocks.size()]);
        double[] xDatadouble = new double[xDataArray.length];
        for (int i = 0; i < xDataArray.length; i++) {
            xDatadouble[i] = xDataArray[i];
        }
        return xDatadouble;
    }

    public List<List<Double>> blockWithIndex(int count, int index, List<List<Double>> dataError) {

        List<Double> xData = new ArrayList<Double>(dataError.get(0));
        List<Double> yData = new ArrayList<Double>(dataError.get(1));

        int start = RANDOM_COUNT / count * index;
        int end = RANDOM_COUNT / count * (index + 1);

        xData = xData.subList(start, end);
        yData = yData.subList(start, end);

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);
        return data;
    }

    public List<List<Double>> blocksWithoutIndex(int count, int index, List<List<Double>> dataError) {

        List<Double> xData = new ArrayList<Double>(dataError.get(0));
        List<Double> yData = new ArrayList<Double>(dataError.get(1));

        int start = RANDOM_COUNT / count * index;
        int end = RANDOM_COUNT / count * (index + 1);

        for (int i = start; i < end; i++) {
            xData.remove(start);
            yData.remove(start);
        }

        List<List<Double>> data = new ArrayList<List<Double>>();
        data.add(xData);
        data.add(yData);
        return data;
    }


    public static void main(String... aArgs) {

        MNK mnk = new MNK();

        List<List<Double>> dataError = mnk.funcDataError;
        List<Double> xData = dataError.get(0);
        List<Double> yData = dataError.get(1);

        List<Chart> charts = new ArrayList<Chart>();
        Chart chart = new Chart(800, 600);
        chart.getStyleManager().setYAxisMin(BOTTOM_Y);
        chart.getStyleManager().setYAxisMax(TOP_Y);
        chart.getStyleManager().setXAxisMax(RIGHT_X);
        chart.getStyleManager().setXAxisMin(LEFT_X);
        chart.getStyleManager().setChartType(ChartType.Scatter);
        chart.getStyleManager().setChartType(ChartType.Scatter);
        chart.getStyleManager().setChartTitleVisible(false);
        chart.getStyleManager().setLegendPosition(LegendPosition.OutsideE);
        chart.getStyleManager().setMarkerSize(10);

        chart.addSeries("Gaussian Blob", xData, yData);
        charts.add(chart); // POINT CHART

        Double[] xDataArray = xData.toArray(new Double[xData.size()]);
        Double[] yDataArray = yData.toArray(new Double[yData.size()]);

        double[] xDatadouble = new double[xDataArray.length];
        double[] yDatadouble = new double[yDataArray.length];

        for (int i = 0; i < xDataArray.length; i++) {
            xDatadouble[i] = xDataArray[i];
        }

        for (int i = 0; i < yDataArray.length; i++) {
            yDatadouble[i] = yDataArray[i];
        }


        double[][] xyTable = {xDatadouble, yDatadouble};
        charts.add(mnk.chartForPower(3, xyTable));
        charts.add(mnk.chartForPower(4, xyTable));
        charts.add(mnk.chartForPower(5, xyTable));
        charts.add(mnk.chartForPower(6, xyTable));
        charts.add(mnk.chartForPower(7, xyTable));
        charts.add(mnk.chartForPower(8, xyTable));
        charts.add(mnk.chartForPower(9, xyTable));
        charts.add(mnk.chartForPower(10, xyTable));

        new SwingWrapper(charts).displayChartMatrix();

        mnk.crossValidation();
    }

    Chart chartForPower(Integer power, double[][] xyTable) {
        double[] result = LeastSquares(xyTable, power);

        List<Double> xDataResult = new ArrayList<Double>();
        List<Double> yDataResult = new ArrayList<Double>();

        for (double x = LEFT_X; x <= RIGHT_X; x += RANGE_X) {
            double y = fy(result, x);
            yDataResult.add(y);
            xDataResult.add(x);
        }

        Chart chart = new Chart(800, 600);
        chart.getStyleManager().setYAxisMin(BOTTOM_Y);
        chart.getStyleManager().setYAxisMax(TOP_Y);
        chart.getStyleManager().setXAxisMax(RIGHT_X);
        chart.getStyleManager().setXAxisMin(LEFT_X);

        chart.addSeries(power.toString(), xDataResult, yDataResult);

        List<List<Double>> data = funcData;
        List<Double> xData = data.get(0);
        List<Double> yData = data.get(1);

        chart.addSeries("V", xData, yData);

        return chart;
    }
}
