using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using FFTClass;



namespace WindowsFormsApplication2
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void dataGridView1_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e){

            int N = dataGridView1.RowCount - 1;


            List<double> xnList = new List<double>(N);
            for (int i = 0; i < N; ++i){

                xnList.Add((Double.Parse(dataGridView1.Rows[i].Cells[0].Value.ToString())));

            }
            double[] xn = xnList.ToArray();



            double Fs = 1;


            FFT fft = new FFT(xn , Fs);
            // List<double> Sol = fft.CalcAmp();



            MessageBox.Show("끝");



            for (int i = 0; i < N / 2 -1 ; ++i){
                // dataGridView1.Rows[i].Cells[1].Value = Sol[i].ToString();
            }


        }



    }
}
