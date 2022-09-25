using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO.Ports;
using System.Threading;
using System.IO;
using System.Globalization;

namespace compass_calibrator
{
    public partial class Form1 : Form
    {
        static string string_from_usart = "";
        static Boolean port_close_flag;
        double X_serial_value, Y_serial_value, Z_serial_value;
        int empty_serial_data_counter;
        string curDir;
        string symbol_mask = "1234567890.-";

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            //change regional standard to en-US
            CultureInfo en = new CultureInfo("en-US");
            Thread.CurrentThread.CurrentCulture = en;

            //current dir name
            curDir = System.IO.Path.GetDirectoryName(
            System.Reflection.Assembly.GetExecutingAssembly().GetModules()[0].FullyQualifiedName);

            stop_button.Enabled = false;
            point_buttons_enable(false);
        }

        private void point_buttons_enable(Boolean b_enable)
        {
            if (b_enable)
            {
                buttonXplus_0.Enabled = true;
                buttonXplus_180.Enabled = true;
                buttonXminus_0.Enabled = true;
                buttonXminus_180.Enabled = true;
                buttonYplus_0.Enabled = true;
                buttonYplus_180.Enabled = true;
                buttonYminus_0.Enabled = true;
                buttonYminus_180.Enabled = true;
                buttonZplus_0.Enabled = true;
                buttonZplus_180.Enabled = true;
                buttonZminus_0.Enabled = true;
                buttonZminus_180.Enabled = true;
            }
            else
            {
                buttonXplus_0.Enabled = false;
                buttonXplus_180.Enabled = false;
                buttonXminus_0.Enabled = false;
                buttonXminus_180.Enabled = false;
                buttonYplus_0.Enabled = false;
                buttonYplus_180.Enabled = false;
                buttonYminus_0.Enabled = false;
                buttonYminus_180.Enabled = false;
                buttonZplus_0.Enabled = false;
                buttonZplus_180.Enabled = false;
                buttonZminus_0.Enabled = false;
                buttonZminus_180.Enabled = false;
            }
        }

        private void start_button_Click(object sender, EventArgs e)
        {
            try
            {
                port.PortName = comboBox1.SelectedItem.ToString();
                //port opening
                port_close_flag = false;
                port.Open();
                stop_button.Enabled = true;
                point_buttons_enable(true);
                start_button.Enabled = false;
                comboBox1.Enabled = false;
                //for the handler
                port.DataReceived += new SerialDataReceivedEventHandler(DataReceivedHandler);
                //timer activation
                timer1.Enabled = true;
            }
            catch
            {
                close_port();
                MessageBox.Show("Empty port name.", "Warning!");
            }
        }

        private void stop_button_Click(object sender, EventArgs e)
        {
            close_port();
        }

        void close_port()
        {
            //timet deactivation
            timer1.Enabled = false;
            //port closing
            port_close_flag = true;
            stop_button.Enabled = false;
            point_buttons_enable(false);
            Thread.Sleep(500);
            port.Close();
            start_button.Enabled = true;
            comboBox1.Enabled = true;
            empty_serial_data_counter = 0;
            string_from_usart = "";
            X_serial_value = 0;
            Y_serial_value = 0;
            Z_serial_value = 0;
            Xlabel.Text = "X = " + X_serial_value.ToString("0.###");
            Ylabel.Text = "Y = " + Y_serial_value.ToString("0.###");
            Zlabel.Text = "Z = " + Z_serial_value.ToString("0.###");
        }

        private static void DataReceivedHandler(object sender, SerialDataReceivedEventArgs e)
        {
            if (!port_close_flag)
            {
                SerialPort sp = (SerialPort)sender;
                string indata = sp.ReadLine();
                string_from_usart = indata;
            }
        }

        private void var_refresh()
        {
            try
            {
                if (string_from_usart != "")
                {
                    string[] serial_values = string_from_usart.Split(',');
                    if (serial_values[0] != "" && serial_values[1] != "" && serial_values[2] != "")
                    {
                        X_serial_value = double.Parse(serial_values[0]);
                        Y_serial_value = double.Parse(serial_values[1]);
                        Z_serial_value = double.Parse(serial_values[2]);
                    }
                }
                else
                {
                    empty_serial_data_counter++;
                    if (empty_serial_data_counter >= 10)
                    {
                        close_port();
                        MessageBox.Show("No serial data.", "Warning!");
                    }
                }
            }
            catch
            {
                close_port();
                MessageBox.Show("Serial port error.", "Warning!");
            }
        }

        private void indication()
        {
            Xlabel.Text = "X = " + X_serial_value.ToString("0.###");
            Ylabel.Text = "Y = " + Y_serial_value.ToString("0.###");
            Zlabel.Text = "Z = " + Z_serial_value.ToString("0.###");
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            var_refresh();
            indication();
        }

        private void comboBox1_DropDown(object sender, EventArgs e)
        {
            comboBox1.Items.Clear();
            foreach (string s in SerialPort.GetPortNames())
            {
                comboBox1.Items.Add(s);
            }    
        }

        private void Form1_FormClosing(object sender, FormClosingEventArgs e)
        {
            close_port();
        }

        private void buttonXplus_0_Click(object sender, EventArgs e)
        {
            textBoxXplus_X_0.Text = X_serial_value.ToString("0.###");
            textBoxXplus_Y_0.Text = Y_serial_value.ToString("0.###");
            textBoxXplus_Z_0.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonXplus_180_Click(object sender, EventArgs e)
        {
            textBoxXplus_X_180.Text = X_serial_value.ToString("0.###");
            textBoxXplus_Y_180.Text = Y_serial_value.ToString("0.###");
            textBoxXplus_Z_180.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonXminus_0_Click(object sender, EventArgs e)
        {
            textBoxXminus_X_0.Text = X_serial_value.ToString("0.###");
            textBoxXminus_Y_0.Text = Y_serial_value.ToString("0.###");
            textBoxXminus_Z_0.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonXminus_180_Click(object sender, EventArgs e)
        {
            textBoxXminus_X_180.Text = X_serial_value.ToString("0.###");
            textBoxXminus_Y_180.Text = Y_serial_value.ToString("0.###");
            textBoxXminus_Z_180.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonYplus_0_Click(object sender, EventArgs e)
        {
            textBoxYplus_X_0.Text = X_serial_value.ToString("0.###");
            textBoxYplus_Y_0.Text = Y_serial_value.ToString("0.###");
            textBoxYplus_Z_0.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonYplus_180_Click(object sender, EventArgs e)
        {
            textBoxYplus_X_180.Text = X_serial_value.ToString("0.###");
            textBoxYplus_Y_180.Text = Y_serial_value.ToString("0.###");
            textBoxYplus_Z_180.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonZplus_0_Click(object sender, EventArgs e)
        {
            textBoxZplus_X_0.Text = X_serial_value.ToString("0.###");
            textBoxZplus_Y_0.Text = Y_serial_value.ToString("0.###");
            textBoxZplus_Z_0.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonZplus_180_Click(object sender, EventArgs e)
        {
            textBoxZplus_X_180.Text = X_serial_value.ToString("0.###");
            textBoxZplus_Y_180.Text = Y_serial_value.ToString("0.###");
            textBoxZplus_Z_180.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonYminus_0_Click(object sender, EventArgs e)
        {
            textBoxYminus_X_0.Text = X_serial_value.ToString("0.###");
            textBoxYminus_Y_0.Text = Y_serial_value.ToString("0.###");
            textBoxYminus_Z_0.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonYminus_180_Click(object sender, EventArgs e)
        {
            textBoxYminus_X_180.Text = X_serial_value.ToString("0.###");
            textBoxYminus_Y_180.Text = Y_serial_value.ToString("0.###");
            textBoxYminus_Z_180.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonZminus_0_Click(object sender, EventArgs e)
        {
            textBoxZminus_X_0.Text = X_serial_value.ToString("0.###");
            textBoxZminus_Y_0.Text = Y_serial_value.ToString("0.###");
            textBoxZminus_Z_0.Text = Z_serial_value.ToString("0.###");
        }

        private void buttonZminus_180_Click(object sender, EventArgs e)
        {
            textBoxZminus_X_180.Text = X_serial_value.ToString("0.###");
            textBoxZminus_Y_180.Text = Y_serial_value.ToString("0.###");
            textBoxZminus_Z_180.Text = Z_serial_value.ToString("0.###");
        }

        private void CalculateButton_Click(object sender, EventArgs e)
        {
            try
            {
                calculate_transformation_matrix();
            }
            catch
            {
                textBox_matrixX_x.Clear();
                textBox_matrixX_y.Clear();
                textBox_matrixX_z.Clear();
                textBox_matrixY_x.Clear();
                textBox_matrixY_y.Clear();
                textBox_matrixY_z.Clear();
                textBox_matrixZ_x.Clear();
                textBox_matrixZ_y.Clear();
                textBox_matrixZ_z.Clear();
                textBox_biasX.Clear();
                textBox_biasY.Clear();
                textBox_biasZ.Clear();
                MessageBox.Show("Wrong input data!", "Warning!");          
            }
        }

        private void calculate_transformation_matrix()
        {
            //Axis X--------------------------------------------------------------------------------------------------
            double[] Xplus_center = new double[3];
            //Centers of the circles
            Xplus_center[0] = (double.Parse(textBoxXplus_X_0.Text) + double.Parse(textBoxXplus_X_180.Text)) / 2;
            Xplus_center[1] = (double.Parse(textBoxXplus_Y_0.Text) + double.Parse(textBoxXplus_Y_180.Text)) / 2;
            Xplus_center[2] = (double.Parse(textBoxXplus_Z_0.Text) + double.Parse(textBoxXplus_Z_180.Text)) / 2;
            //Centers of the circles
            double[] Xminus_center = new double[3];
            Xminus_center[0] = (double.Parse(textBoxXminus_X_0.Text) + double.Parse(textBoxXminus_X_180.Text)) / 2;
            Xminus_center[1] = (double.Parse(textBoxXminus_Y_0.Text) + double.Parse(textBoxXminus_Y_180.Text)) / 2;
            Xminus_center[2] = (double.Parse(textBoxXminus_Z_0.Text) + double.Parse(textBoxXminus_Z_180.Text)) / 2;
            //Vector from the center of minus circle to the center of plus circle
            double[] Xvector = new double[3];
            Xvector[0] = Xplus_center[0] - Xminus_center[0];
            Xvector[1] = Xplus_center[1] - Xminus_center[1];
            Xvector[2] = Xplus_center[2] - Xminus_center[2];

            //Axis Y--------------------------------------------------------------------------------------------------
            double[] Yplus_center = new double[3];
            //Centers of the circles
            Yplus_center[0] = (double.Parse(textBoxYplus_X_0.Text) + double.Parse(textBoxYplus_X_180.Text)) / 2;
            Yplus_center[1] = (double.Parse(textBoxYplus_Y_0.Text) + double.Parse(textBoxYplus_Y_180.Text)) / 2;
            Yplus_center[2] = (double.Parse(textBoxYplus_Z_0.Text) + double.Parse(textBoxYplus_Z_180.Text)) / 2;
            //Centers of the circles
            double[] Yminus_center = new double[3];
            Yminus_center[0] = (double.Parse(textBoxYminus_X_0.Text) + double.Parse(textBoxYminus_X_180.Text)) / 2;
            Yminus_center[1] = (double.Parse(textBoxYminus_Y_0.Text) + double.Parse(textBoxYminus_Y_180.Text)) / 2;
            Yminus_center[2] = (double.Parse(textBoxYminus_Z_0.Text) + double.Parse(textBoxYminus_Z_180.Text)) / 2;
            //Vector from the center of minus circle to the center of plus circle
            double[] Yvector = new double[3];
            Yvector[0] = Yplus_center[0] - Yminus_center[0];
            Yvector[1] = Yplus_center[1] - Yminus_center[1];
            Yvector[2] = Yplus_center[2] - Yminus_center[2];

            //Axis Z--------------------------------------------------------------------------------------------------
            double[] Zplus_center = new double[3];
            //Centers of the circles
            Zplus_center[0] = (double.Parse(textBoxZplus_X_0.Text) + double.Parse(textBoxZplus_X_180.Text)) / 2;
            Zplus_center[1] = (double.Parse(textBoxZplus_Y_0.Text) + double.Parse(textBoxZplus_Y_180.Text)) / 2;
            Zplus_center[2] = (double.Parse(textBoxZplus_Z_0.Text) + double.Parse(textBoxZplus_Z_180.Text)) / 2;
            //Centers of the circles
            double[] Zminus_center = new double[3];
            Zminus_center[0] = (double.Parse(textBoxZminus_X_0.Text) + double.Parse(textBoxZminus_X_180.Text)) / 2;
            Zminus_center[1] = (double.Parse(textBoxZminus_Y_0.Text) + double.Parse(textBoxZminus_Y_180.Text)) / 2;
            Zminus_center[2] = (double.Parse(textBoxZminus_Z_0.Text) + double.Parse(textBoxZminus_Z_180.Text)) / 2;
            //Vector from the center of minus circle to the center of plus circle
            double[] Zvector = new double[3];
            Zvector[0] = Zplus_center[0] - Zminus_center[0];
            Zvector[1] = Zplus_center[1] - Zminus_center[1];
            Zvector[2] = Zplus_center[2] - Zminus_center[2];

            // Rotation matrix--------------------------------------------------------------------------------------
            // rotation_matrix[a][b], a - number of the rows, b - number of the columbs
            double[][] rotation_matrix = new double[3][];
            rotation_matrix[0] = new double[3];
            rotation_matrix[1] = new double[3];
            rotation_matrix[2] = new double[3];
            //Deviding by main value, for example for X axis - deviding by X coordinate, for Y axis by Y coordinate, for Z axis by Z cordinate
            rotation_matrix[0][0] = Xvector[0] / Xvector[0]; rotation_matrix[0][1] = Yvector[0] / Yvector[1]; rotation_matrix[0][2] = Zvector[0] / Zvector[2];
            rotation_matrix[1][0] = Xvector[1] / Xvector[0]; rotation_matrix[1][1] = Yvector[1] / Yvector[1]; rotation_matrix[1][2] = Zvector[1] / Zvector[2];
            rotation_matrix[2][0] = Xvector[2] / Xvector[0]; rotation_matrix[2][1] = Yvector[2] / Yvector[1]; rotation_matrix[2][2] = Zvector[2] / Zvector[2];
            //Matrix inversion
            rotation_matrix = InvertMatrix(rotation_matrix);

            //Determinating of the corrected by ratation matrix centers of the circles 
            Xplus_center = MatrixVectorMultiply(rotation_matrix, Xplus_center);
            Xminus_center = MatrixVectorMultiply(rotation_matrix, Xminus_center);
            Yplus_center = MatrixVectorMultiply(rotation_matrix, Yplus_center);
            Yminus_center = MatrixVectorMultiply(rotation_matrix, Yminus_center);
            Zplus_center = MatrixVectorMultiply(rotation_matrix, Zplus_center);
            Zminus_center = MatrixVectorMultiply(rotation_matrix, Zminus_center);

            //Determinating of the elipsoid center---------------------------------------------------------------------------
            double[] center = new double[3];
            center[0] = (Xplus_center[0] + Xminus_center[0] + Yplus_center[0] + Yminus_center[0] + Zplus_center[0] + Zminus_center[0]) / 6;
            center[1] = (Xplus_center[1] + Xminus_center[1] + Yplus_center[1] + Yminus_center[1] + Zplus_center[1] + Zminus_center[1]) / 6;
            center[2] = (Xplus_center[2] + Xminus_center[2] + Yplus_center[2] + Yminus_center[2] + Zplus_center[2] + Zminus_center[2]) / 6;

            //Determinating of the radius of the future sphere-----------------------------------------------------------------------
            double x_length = Math.Abs(Xplus_center[0] - Xminus_center[0])/2;
            double y_length = Math.Abs(Yplus_center[1] - Yminus_center[1])/2;
            double z_length = Math.Abs(Zplus_center[2] - Zminus_center[2])/2;
            double[] Xplus_0 = new double[3];
            Xplus_0[0] = double.Parse(textBoxXplus_X_0.Text); 
            Xplus_0[1] = double.Parse(textBoxXplus_Y_0.Text); 
            Xplus_0[2] = double.Parse(textBoxXplus_Z_0.Text);
            Xplus_0 = MatrixVectorMultiply(rotation_matrix, Xplus_0);
            double[] Yplus_0 = new double[3];
            Yplus_0[0] = double.Parse(textBoxYplus_X_0.Text);
            Yplus_0[1] = double.Parse(textBoxYplus_Y_0.Text);
            Yplus_0[2] = double.Parse(textBoxYplus_Z_0.Text);
            Yplus_0 = MatrixVectorMultiply(rotation_matrix, Yplus_0);
            double[] Zplus_0 = new double[3];
            Zplus_0[0] = double.Parse(textBoxZplus_X_0.Text);
            Zplus_0[1] = double.Parse(textBoxZplus_Y_0.Text);
            Zplus_0[2] = double.Parse(textBoxZplus_Z_0.Text);
            Zplus_0 = MatrixVectorMultiply(rotation_matrix, Zplus_0);
            double x_abs = Math.Sqrt(x_length * x_length + Xplus_0[1] * Xplus_0[1] + Xplus_0[2] * Xplus_0[2]);
            double y_abs = Math.Sqrt(Yplus_0[0] * Yplus_0[0] + y_length * y_length + Yplus_0[2] * Yplus_0[2]);
            double z_abs = Math.Sqrt(Zplus_0[0] * Zplus_0[0] + Zplus_0[1] * Zplus_0[1] + z_length * z_length);
            //sphere radius
            double sphere_radius = (x_abs + y_abs + z_abs) / 3;

            //Scales for the each axis------------------------------------------------
            //Diameter of the sphere
            double diameter = sphere_radius * 2;
            double kx = Math.Abs(diameter / (Xplus_center[0] - Xminus_center[0]));
            double ky = Math.Abs(diameter / (Yplus_center[1] - Yminus_center[1]));
            double kz = Math.Abs(diameter / (Zplus_center[2] - Zminus_center[2]));

            //Multiplying elements of matrix by scales
            rotation_matrix[0][0] = rotation_matrix[0][0] * kx; rotation_matrix[0][1] = rotation_matrix[0][1] * ky; rotation_matrix[0][2] = rotation_matrix[0][2] * kz;
            rotation_matrix[1][0] = rotation_matrix[1][0] * kx; rotation_matrix[1][1] = rotation_matrix[1][1] * ky; rotation_matrix[1][2] = rotation_matrix[1][2] * kz;
            rotation_matrix[2][0] = rotation_matrix[2][0] * kx; rotation_matrix[2][1] = rotation_matrix[2][1] * ky; rotation_matrix[2][2] = rotation_matrix[2][2] * kz;
            
            //Bias
            double[] bias = new double[3];
            bias[0] = center[0];
            bias[1] = center[1];
            bias[2] = center[2];

            //Indication
            //Transformation matrix
            textBox_matrixX_x.Text = rotation_matrix[0][0].ToString("0.###"); textBox_matrixY_x.Text = rotation_matrix[0][1].ToString("0.###"); textBox_matrixZ_x.Text = rotation_matrix[0][2].ToString("0.###");
            textBox_matrixX_y.Text = rotation_matrix[1][0].ToString("0.###"); textBox_matrixY_y.Text = rotation_matrix[1][1].ToString("0.###"); textBox_matrixZ_y.Text = rotation_matrix[1][2].ToString("0.###");
            textBox_matrixX_z.Text = rotation_matrix[2][0].ToString("0.###"); textBox_matrixY_z.Text = rotation_matrix[2][1].ToString("0.###"); textBox_matrixZ_z.Text = rotation_matrix[2][2].ToString("0.###");
            //Bias
            textBox_biasX.Text = bias[0].ToString("0.###");
            textBox_biasY.Text = bias[1].ToString("0.###");
            textBox_biasZ.Text = bias[2].ToString("0.###");
        }

        public static double[] MatrixVectorMultiply(double[][] matrixA, double[] vectorB)
        {
            int aRows = matrixA.Length; int aCols = matrixA[0].Length;
            int bRows = vectorB.Length;
            if (aCols != bRows)
                throw new Exception("Non-conformable matrices in MatrixProduct");
            double[] result = new double[aRows];
            for (int i = 0; i < aRows; ++i) // each row of A
                for (int k = 0; k < aCols; ++k)
                    result[i] += matrixA[i][k] * vectorB[k];
            return result;
        }

        public static double[][] InvertMatrix(double[][] A)
        {
            int n = A.Length;
            //e will represent each column in the identity matrix
            double[] e;
            //x will hold the inverse matrix to be returned
            double[][] x = new double[n][];
            for (int i = 0; i < n; i++)
            {
                x[i] = new double[A[i].Length];
            }
            /*
            * solve will contain the vector solution for the LUP decomposition as we solve
            * for each vector of x.  We will combine the solutions into the double[][] array x.
            * */
            double[] solve;

            //Get the LU matrix and P matrix (as an array)
            Tuple<double[][], int[]> results = LUPDecomposition(A);

            double[][] LU = results.Item1;
            int[] P = results.Item2;

            /*
            * Solve AX = e for each column ei of the identity matrix using LUP decomposition
            * */
            for (int i = 0; i < n; i++)
            {
                e = new double[A[i].Length];
                e[i] = 1;
                solve = LUPSolve(LU, P, e);
                for (int j = 0; j < solve.Length; j++)
                {
                    x[j][i] = solve[j];
                }
            }
            return x;
        }

        public static double[] LUPSolve(double[][] LU, int[] pi, double[] b)
        {
            int n = LU.Length - 1;
            double[] x = new double[n + 1];
            double[] y = new double[n + 1];
            double suml = 0;
            double sumu = 0;
            double lij = 0;

            /*
            * Solve for y using formward substitution
            * */
            for (int i = 0; i <= n; i++)
            {
                suml = 0;
                for (int j = 0; j <= i - 1; j++)
                {
                    /*
                    * Since we've taken L and U as a singular matrix as an input
                    * the value for L at index i and j will be 1 when i equals j, not LU[i][j], since
                    * the diagonal values are all 1 for L.
                    * */
                    if (i == j)
                    {
                        lij = 1;
                    }
                    else
                    {
                        lij = LU[i][j];
                    }
                    suml = suml + (lij * y[j]);
                }
                y[i] = b[pi[i]] - suml;
            }
            //Solve for x by using back substitution
            for (int i = n; i >= 0; i--)
            {
                sumu = 0;
                for (int j = i + 1; j <= n; j++)
                {
                    sumu = sumu + (LU[i][j] * x[j]);
                }
                x[i] = (y[i] - sumu) / LU[i][i];
            }
            return x;
        }

        public static Tuple<double[][], int[]> LUPDecomposition(double[][] A)
        {
            int n = A.Length - 1;
            /*
            * pi represents the permutation matrix.  We implement it as an array
            * whose value indicates which column the 1 would appear.  We use it to avoid 
            * dividing by zero or small numbers.
            * */
            int[] pi = new int[n + 1];
            double p = 0;
            int kp = 0;
            int pik = 0;
            int pikp = 0;
            double aki = 0;
            double akpi = 0;

            //Initialize the permutation matrix, will be the identity matrix
            for (int j = 0; j <= n; j++)
            {
                pi[j] = j;
            }

            for (int k = 0; k <= n; k++)
            {
                /*
                * In finding the permutation matrix p that avoids dividing by zero
                * we take a slightly different approach.  For numerical stability
                * We find the element with the largest 
                * absolute value of those in the current first column (column k).  If all elements in
                * the current first column are zero then the matrix is singluar and throw an
                * error.
                * */
                p = 0;
                for (int i = k; i <= n; i++)
                {
                    if (Math.Abs(A[i][k]) > p)
                    {
                        p = Math.Abs(A[i][k]);
                        kp = i;
                    }
                }
                if (p == 0)
                {
                    throw new Exception("singular matrix");
                }
                /*
                * These lines update the pivot array (which represents the pivot matrix)
                * by exchanging pi[k] and pi[kp].
                * */
                pik = pi[k];
                pikp = pi[kp];
                pi[k] = pikp;
                pi[kp] = pik;

                /*
                * Exchange rows k and kpi as determined by the pivot
                * */
                for (int i = 0; i <= n; i++)
                {
                    aki = A[k][i];
                    akpi = A[kp][i];
                    A[k][i] = akpi;
                    A[kp][i] = aki;
                }

                /*
                    * Compute the Schur complement
                    * */
                for (int i = k + 1; i <= n; i++)
                {
                    A[i][k] = A[i][k] / A[k][k];
                    for (int j = k + 1; j <= n; j++)
                    {
                        A[i][j] = A[i][j] - (A[i][k] * A[k][j]);
                    }
                }
            }
            return Tuple.Create(A, pi);
        }

        private void help_image_form(string help_image_massge, string help_image_link)
        {
            axis_image_view frm_op = new axis_image_view();
            frm_op.Owner = this;
            frm_op.text = help_image_massge;
            frm_op.image_link = help_image_link;
            frm_op.ShowDialog();
        }

        private void button9_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 0°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "X_plus_0.png");
        }

        private void button10_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 180°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "X_plus_180.png");
        }

        private void button11_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 0°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Y_plus_0.png");
        }

        private void button12_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 180°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Y_plus_180.png");
        }

        private void button13_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 0°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Z_plus_0.png");
        }

        private void button14_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 180°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Z_plus_180.png");
        }

        private void button3_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 0°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "X_minus_0.png");
        }

        private void button4_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 180°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "X_minus_180.png");
        }

        private void button5_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 0°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Y_minus_0.png");
        }

        private void button6_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 180°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Y_minus_180.png");
        }

        private void button7_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 0°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Z_minus_0.png");
        }

        private void button8_Click(object sender, EventArgs e)
        {
            help_image_form("and click \"Point 180°\" button or enter data manually.", curDir + "\\MagMaster Files\\images\\" + "Z_minus_180.png");
        }

        private void help_button_Click(object sender, EventArgs e)
        {
            serial_port_help_from frm_op = new serial_port_help_from();
            frm_op.Owner = this;
            frm_op.help_text = System.IO.File.ReadAllText(curDir + "\\MagMaster Files\\texts\\" + "sphelp.txt"); 
            frm_op.ShowDialog();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            help_how_to_use frm_op = new help_how_to_use();
            frm_op.Owner = this;
            frm_op.image_link = curDir + "\\MagMaster Files\\images\\" + "matrix.png";
            frm_op.ShowDialog();
        }

        private void textBoxXplus_X_0_KeyPress(object sender, KeyPressEventArgs e)
        {
            TextBox tx = (TextBox)sender;
            if ((symbol_mask.IndexOf(e.KeyChar) != -1) || (e.KeyChar == 8))
            {
                if ((e.KeyChar == '.') && (tx.Text.IndexOf(".") != -1))
                    e.Handled = true;
                if ((e.KeyChar == '-') && (tx.Text.IndexOf("-") != -1))
                    e.Handled = true;
            }
            else
                e.Handled = true;       
        }

        private void textBoxXplus_X_0_TextChanged(object sender, EventArgs e)
        {
            TextBox tx = (TextBox)sender;
            string strbox = tx.Text;

            string symbol_st = ".";
            if (strbox.Length >= symbol_st.Length)
            {
                string substrbox = strbox.Substring(0, symbol_st.Length);
                if (substrbox == symbol_st) strbox = "" + strbox.Substring(symbol_st.Length, strbox.Length - symbol_st.Length);
            }

            symbol_st = "-.";
            if (strbox.Length >= symbol_st.Length)
            {
                string substrbox = strbox.Substring(0, symbol_st.Length);
                if (substrbox == symbol_st) strbox = "-" + strbox.Substring(symbol_st.Length, strbox.Length - symbol_st.Length);
            }

            symbol_st = "-";
            if (strbox.Length >= symbol_st.Length)
            {
                string substrbox = strbox.Substring(0, symbol_st.Length);
                if (substrbox != symbol_st) strbox = strbox.Replace("-", "");
            }

            tx.Text = strbox;

        }

        private void button2_Click(object sender, EventArgs e)
        {
            System.Diagnostics.Process.Start("http://diydrones.com/profiles/blogs/advanced-hard-and-soft-iron-magnetometer-calibration-for-dummies");
        }
    }
}
