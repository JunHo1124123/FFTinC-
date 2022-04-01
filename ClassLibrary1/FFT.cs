using System;
using System.Collections.Generic;
using System.Collections;
using System.Numerics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FFTClass{
    public class FFT{

        private Complex[] cpx_xn;
        private double[] xn;
        private int N;
        private double Fs;
        public FFT(double[] xn_, double Fs_){

            
            xn = xn_;
            N = xn_.Length;
            Fs = Fs_;



            // IntegerFactorization(N);


            List<Complex> cpx_xn_List = new List<Complex>(N);
            for (int i = 0; i < N; ++i){

                cpx_xn_List.Add(new Complex(xn[i], 0));
            }
            cpx_xn = cpx_xn_List.ToArray();

        }







        private int n_stages;    // 계산할 radix 횟수 = 반복 계산 횟수 - 1
        private List<int> stage;    // redix order 백터
        private void IntegerFactorization(int num)
        {


            // 한번에 연산 할 최대 radix 값
            // 반복적으로 확인이 필요함.
            // radix 2^n 개 기준으로, n=4까지는 FFT 계산보다 DFT 계산이 빠름 (논문 참고)
            // 하지만 알고리즘 및 하드웨어의 차이가 있어, 실험적으로 접근이 필요함.
            // 안 넣는 것 보다는 20~30사이의 값을 입력하는 것이 빠른 것을 확인함
            int Max_radix = 24;



            // 소인수 분해 알고리즘
            // 작은 숫자부터 큰 숫자 순으로 분해
            // 1은 생략함.
            List<int> origin_stage = new List<int>();
            int factor = 2;
            while (num != 1)
            {
                if (num % factor == 0)
                {
                    origin_stage.Add(factor);
                    num /= factor;
                }
                else
                {
                    factor++;
                }
            }
            int n_origin_stage = origin_stage.Count();



            // 소인수 분해를 적절한 사이즈의 radix로 계산하기 위해 합치는 알고리즘
            // 최대 radix order가 Max_radix가 되도록 계산
            // 끝값과 처음값들을 곱함.
            // 그 값이 Max_radix보다 작으면 곱한 상태로 stage 배열에 입력, 크면 곱하지 않고 stage 배열에 입력
            // stage 개수를 줄여서 연산 속도 증가
            int Factor_mul;
            int remove_j = 0;

            for (int i = n_origin_stage - 1; i - remove_j >= 0; --i)
            {
                Factor_mul = origin_stage[i];           // *(origin_stage_ptr + i);
                for (int j = 0 + remove_j; j <= i; ++j)
                {
                    Factor_mul *= origin_stage[j];                 // *(origin_stage_ptr + j);
                    if (Factor_mul > Max_radix || i == j)
                    {
                        remove_j = j;
                        stage.Add(Factor_mul / origin_stage[j]);  // *(origin_stage_ptr + j);
                        break;
                    }
                }
            }
            n_stages = stage.Count();



        }













        private List<int> Mapping_idx;
        private void MappingIndexing(int Nx, int Ny)
        {

            Mapping_idx.Capacity = N;

            int Ni = Nx * Ny;
            int ky, kx, n;

            for (int k = 0; k < N; ++k)
            {

                ky = k % Ny;
                kx = k / Ny;

                Mapping_idx[k] = (kx % Nx) + (kx / Nx) * Ni + ky * Nx;

                // *(Mapping_idx_ptr + k) = (kx % Nx) + (kx / Nx) * Ni + ky * Nx;

            }


        }




        private List<int> digit_reverse_idx;
        private void Digit_Reverse(List<int> stage)
        {

            digit_reverse_idx.Capacity = N;

            int k, Nx, Ny, Ni;

            // 모든 데이터 n개에 대해서 모두 수행
            for (int n = 0; n < N; ++n)
            {

                // Stage 0에서의 초기 위치 및 초기 Nx
                k = n;
                Nx = stage[0];  // *(stage_ptr);

                // Stage 1부터 끝까지 인덱스 계산 --> 마지막 인덱스 확인 가능
                for (int s = 1; s < n_stages; ++s)
                {

                    Ny = stage[s];  // *(stage_ptr + s);
                    Ni = Ny * Nx;

                    // s 번째 스테이지에서의 인덱스 업데이트
                    k = (k * Ny) % Ni + (k / Nx) % Ny + Ni * (k / Ni);

                    // s+1 번째 스테이지의 계산을 위한 Nx 업데이트
                    Nx *= Ny;
                }

                // 계산 완료된 인덱스값을 배열에 입력

                digit_reverse_idx[n] = k;
                //*(digit_reverse_idx_ptr + n) = k;

            }




        }



        private List<Complex> make_W_vec(int Nx, int Ny)
        {

            List<Complex> W_vec = new List<Complex>();
            W_vec.Capacity = (Nx - 1) * (Ny - 1) + 1;

            int Ni = Nx * Ny;
            int min = (Nx + Ny - Math.Abs(Nx - Ny)) / 2;
            int max = Nx + Ny - min;

            // 해당 되지 않는 값에 1 입력
            for (int i = 0; i < N; ++i)
            {
                W_vec[i] = new Complex(1, 0);
            }


            // 해당 값에 W값 입력
            for (int nx = 1; nx < min; ++nx)
            {
                for (int ky = nx; ky < max; ++ky)
                {
                    W_vec[ky * nx] = new Complex(
                        Math.Cos(2 * Math.PI * nx * ky / Ni),
                        -1 * Math.Sin(2 * Math.PI * nx * ky / Ni)
                        );

                }
            }

            return W_vec;
        }





        private List<Complex> Solve_DFT;
        private void DFT(int radix, List<Complex> xn)
        {

            Solve_DFT.Clear();
            Solve_DFT.Capacity = radix;
            Complex z;

            // 최적화 필요
            for (int k = 0; k < radix; k++)
            {
                for (int n = 0; n < radix; n++)
                {
                    z = new Complex(
                        Math.Cos(2 * Math.PI * n * k / radix),
                        -1 * Math.Sin(2 * Math.PI * n * k / radix)
                        );
                    Solve_DFT[k] += xn[n] * z;
                    // *(Solve_DFT_1st_pointer + k) += *(xn_1st_pointer + n) * cpx(cos(2 * PI * n * k / radix), -1 * sin(2 * PI * n * k / radix));
                }
            }





        }


        private List<Complex> FFT_Result;
        private void Transform()
        {

            // 데이터 개수 N이 소수일 경우 -> radix-N 계산 (DFT 변환)
            if (n_stages == 1) { DFT(stage[0], cpx_xn); FFT_Result = Solve_DFT; }

            // 데이터 개수 N이 합성수일 경우 -> mixed-radix FFT 수행
            else if (n_stages > 1)
            {

                ///////////////////////////////////////// 소인수 분해 결과로 Digit_Reverse 작업 수행

                List<Complex> Digit_Reversed_Data_vector = null;
                Digit_Reversed_Data_vector.Capacity = N;
                Digit_Reverse(stage);

                for (int i = 0; i < N; ++i)
                {
                    Digit_Reversed_Data_vector[i] = cpx_xn[digit_reverse_idx[i]];
                }

                /////////////////////////////////////////////////////////////////// 첫번째 radix 수행

                List<Complex> Done_DFT_Data = new List<Complex>();
                Done_DFT_Data.Capacity = N;

                List<Complex> Do_DFT_Data = new List<Complex>();
                Do_DFT_Data.Capacity = stage[0];

                for (int i = 0; i < (N / stage[0]); ++i)
                {
                    for (int j = 0; j < stage[0]; ++j)
                    {
                        Do_DFT_Data[j] = Digit_Reversed_Data_vector[i * stage[0] + j];
                    }

                    DFT(stage[0], Do_DFT_Data);

                    for (int j = 0; j < stage[0]; ++j)
                    {
                        Done_DFT_Data[i * stage[0] + j] = Solve_DFT[j];
                    }

                }


                /////////////////////////////////////////////////////////////// radix 수행 루프 생성

                int Nx = stage[0];
                int Ny;

                List<Complex> Reversed_Data = null;
                Reversed_Data.Capacity = N;

                for (int s = 1; s < n_stages; ++s)
                {

                    Ny = stage[s];

                    MappingIndexing(Nx, Ny);


                    //////////////////////////////////////////////// twidle factor 곱 & Reverse 작업

                    List<Complex> W_vec = make_W_vec(Nx, Ny);

                    for (int n = 0; n < N; ++n)
                    {
                        // W_vec 매개변수: (nx*ky)  -->  nx, ky 구하는 공식을 직접 입력함
                        Reversed_Data[n] = Done_DFT_Data[Mapping_idx[n]] * W_vec[(Mapping_idx[n] % Nx) * ((Mapping_idx[n] / Nx) % Ny)];
                    }

                    ///////////////////////////////////////////////////// s번째 radix 수행

                    Done_DFT_Data.Clear();
                    Do_DFT_Data.Clear();
                    Do_DFT_Data.Capacity = stage[0];


                    for (int i = 0; i < (N / stage[s]); ++i)
                    {
                        for (int j = 0; j < stage[s]; ++j)
                        {
                            Do_DFT_Data[j] = Reversed_Data[i * stage[s] + j];
                        }

                        DFT(stage[s], Do_DFT_Data);

                        for (int j = 0; j < stage[s]; ++j)
                        {
                            Done_DFT_Data[i * stage[s] + j] = Solve_DFT[j];
                        }

                    }

                    //////////////////////////////////////////////////////데이터 원복하기

                    for (int n = 0; n < N; ++n)
                    {

                        Reversed_Data[Mapping_idx[n]] = Done_DFT_Data[n];
                        // *(Reversed_ptr + *(Mapping_idx_ptr + n)) = *(Done_DFT_ptr + n);

                    }

                    Nx *= Ny;

                    Done_DFT_Data = Reversed_Data;



                }

                FFT_Result = Done_DFT_Data;

            }





        }

        public List<double> CalcAmp(string window = "None"){

            List<double> Amplitude = new List<double>();
            Amplitude.Capacity = N / 2 + 1;

	        if (window == "None"){}
	        else if (window == "blackman"){

		        for (int i = 0; i < N; ++i){
                    cpx_xn[i] *= 0.42 - (0.5 * Math.Cos(2 * Math.PI * i / N)) + (0.08 * Math.Cos(4 * Math.PI * i / N));
			        // *(cpx_xn_ptr + i) *= 0.42 - (0.5*cos(2 * PI*i / N)) + (0.08*cos(4 * PI*i / N));
		        }

	        }
	        else if (window == "hanning" || window == "hann"){

		        for (int i = 0; i < N; ++i){
                    cpx_xn[i] *= 0.5 * (1 - Math.Cos(2 * Math.PI * i / N));
			        // *(cpx_xn_ptr + i) *= 0.5*(1 - cos(2 * PI*i / N));
		        }

	        }
	        else if (window == "hamming"){

		        for (int i = 0; i < N; ++i){
                    cpx_xn[i] *= 0.54 - (0.46 * Math.Cos(2 * Math.PI * i / N));
			        // *(cpx_xn_ptr + i) *= 0.54 - (0.46*cos(2 * PI*i / N));
		        }

	        }
	        else if (window == "flattop"){

                double a = 2.0 * Math.PI / N;
		        for (int i = 0; i < N; ++i){
                    cpx_xn[i] *= (1) - (1.93) * Math.Cos(a * i) + (1.29) * Math.Cos(2 * a * i) - (0.388) * Math.Cos(3 * a * i) + (0.028) * Math.Cos(4 * a * i);
			        // *(cpx_xn_ptr + i) *= (1) - (1.93)*cos(a*i) + (1.29)*cos(2 * a*i) - (0.388)*cos(3 * a*i) + (0.028)*cos(4 * a*i);
		        }

	        }
	        else if (window == "gaussian"){
		        // 추가 필요

		        for (int i = 0; i < N; ++i){
                    cpx_xn[i] *= 0;
		        }

	        }
	        else{

		        for (int i = 0; i < N; ++i){
                    cpx_xn[i] *= 0;
		        }
		        // 다시 입력할 것
	        }


	        Transform();
	

	        for (int i = 0; i <= N / 2; ++i){
                Amplitude[i] = Complex.Abs(FFT_Result[i]) / N * 2;
	        }


            return Amplitude;
        }

    }
}
