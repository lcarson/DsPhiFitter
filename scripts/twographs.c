void twographs() {
   TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);
   const Int_t n = 20;
   Double_t x[n], y1[n], y2[n];
   for (Int_t i=0;i<n;i++) {
     x[i]  = i*0.1;
     y1[i] = 10*sin(x[i]+0.2);
     y2[i] =  9*sin(x[i]+0.2);
   }
   gr1 = new TGraph(n,x,y1);
   gr1->SetFillColor(2);
   gr2 = new TGraph(n,x,y2);
   gr2->SetFillColor(0);
   gr1->Draw("AF");
   gr2->Draw("F");
}

void grshade() {
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

   c1->SetGrid();
   c1->DrawFrame(0,0,2.2,12);
   
   const Int_t n = 20;
   Double_t x[n], y[n],ymin[n], ymax[n];
   Int_t i;
   for (i=0;i<n;i++) {
     x[i] = 0.1+i*0.1;
     ymax[i] = 10*sin(x[i]+0.2);
     ymin[i] = 8*sin(x[i]+0.1);
     y[i] = 9*sin(x[i]+0.15);
   }
   TGraph *grmin = new TGraph(n,x,ymin);
   TGraph *grmax = new TGraph(n,x,ymax);
   TGraph *gr    = new TGraph(n,x,y);
   TGraph *grshade = new TGraph(2*n);
   for (i=0;i<n;i++) {
      grshade->SetPoint(i,x[i],ymax[i]);
      grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
   }
   grshade->SetFillStyle(3013);
   grshade->SetFillColor(16);
   grshade->Draw("f");
   grmin->Draw("l");
   grmax->Draw("l");
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("CP");
}