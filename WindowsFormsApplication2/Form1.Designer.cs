namespace WindowsFormsApplication2
{
    partial class Form1
    {
        /// <summary>
        /// 필수 디자이너 변수입니다.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// 사용 중인 모든 리소스를 정리합니다.
        /// </summary>
        /// <param name="disposing">관리되는 리소스를 삭제해야 하면 true이고, 그렇지 않으면 false입니다.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form 디자이너에서 생성한 코드

        /// <summary>
        /// 디자이너 지원에 필요한 메서드입니다.
        /// 이 메서드의 내용을 코드 편집기로 수정하지 마십시오.
        /// </summary>
        private void InitializeComponent()
        {
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.BeforeFFT = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.AfterFFT = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.FFTStartBtn = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            this.SuspendLayout();
            // 
            // dataGridView1
            // 
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.BeforeFFT,
            this.AfterFFT});
            this.dataGridView1.Location = new System.Drawing.Point(827, 57);
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.RowTemplate.Height = 27;
            this.dataGridView1.Size = new System.Drawing.Size(351, 535);
            this.dataGridView1.TabIndex = 0;
            this.dataGridView1.CellContentClick += new System.Windows.Forms.DataGridViewCellEventHandler(this.dataGridView1_CellContentClick);
            // 
            // BeforeFFT
            // 
            this.BeforeFFT.HeaderText = "Before FFT";
            this.BeforeFFT.Name = "BeforeFFT";
            this.BeforeFFT.Resizable = System.Windows.Forms.DataGridViewTriState.True;
            this.BeforeFFT.SortMode = System.Windows.Forms.DataGridViewColumnSortMode.NotSortable;
            this.BeforeFFT.Width = 150;
            // 
            // AfterFFT
            // 
            this.AfterFFT.HeaderText = "After FFT";
            this.AfterFFT.Name = "AfterFFT";
            this.AfterFFT.ReadOnly = true;
            this.AfterFFT.Width = 150;
            // 
            // FFTStartBtn
            // 
            this.FFTStartBtn.Location = new System.Drawing.Point(525, 244);
            this.FFTStartBtn.Name = "FFTStartBtn";
            this.FFTStartBtn.Size = new System.Drawing.Size(165, 65);
            this.FFTStartBtn.TabIndex = 1;
            this.FFTStartBtn.Text = "FFT 수행";
            this.FFTStartBtn.UseVisualStyleBackColor = true;
            this.FFTStartBtn.Click += new System.EventHandler(this.button1_Click);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Font = new System.Drawing.Font("굴림", 20F);
            this.label1.Location = new System.Drawing.Point(546, 161);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(110, 34);
            this.label1.TabIndex = 2;
            this.label1.Text = "label1";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 15F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1314, 768);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.FFTStartBtn);
            this.Controls.Add(this.dataGridView1);
            this.Name = "Form1";
            this.Text = "Form1";
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.Button FFTStartBtn;
        private System.Windows.Forms.DataGridViewTextBoxColumn BeforeFFT;
        private System.Windows.Forms.DataGridViewTextBoxColumn AfterFFT;
        private System.Windows.Forms.Label label1;


    }
}

