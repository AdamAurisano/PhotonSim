void vis()
{
  gStyle->SetCanvasPreferGL(true);
  TCanvas* c1 = new TCanvas("c1");
  TGeoManager* manager = new TGeoManager("world", "world");
  TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
  TGeoMedium   *med = new TGeoMedium("Vacuum", 1, mat);

  double D = 0.07;
  double claddingPercent = 0.03;
  double rCore = D*(0.5 - 2*claddingPercent);
  double rInner = rCore + claddingPercent*D;
  double rOuter = rInner + claddingPercent*D;
  
  TGeoVolume *top = manager->MakeTube("Top", med, 0.0, 40, 500);
  manager->SetTopVolume(top);
  TGeoVolume *core  = manager->MakeTube("Core",  med, 0.0*1000,    rCore*1000, 500);
  TGeoVolume *inner = manager->MakeTube("Inner", med, rCore*1000,  rInner*1000, 500);
  TGeoVolume *outer = manager->MakeTube("Outer", med, rInner*1000, rOuter*1000, 500);

  top->AddNode(core,1);
  top->AddNode(inner,1);
  top->AddNode(outer,1);
  
  manager->CloseGeometry();

  core->SetLineColor(kGreen);
  inner->SetLineColor(kBlue);
  outer->SetLineColor(kRed);

  core->SetTransparency(75);
  inner->SetTransparency(75);
  outer->SetTransparency(75);

  
  //top->SetLineColor(kMagenta);
  //top->SetTransparency(50);
  //manager->SetTopVisible();
  top->Draw("ogl");

  string input;
  bool done = false;
  
  TRandom gen;
  TPolyLine3D* poly = new TPolyLine3D();
  poly->SetLineWidth(5);
  poly->SetLineColor(kRed);
  TTimer *timer = new TTimer("gSystem->ProcessEvents();",50,false);

  do {
    timer->TurnOn();
    timer->Reset();
    cout << "Type <return> for next event, <q> to exit: " << endl;
    getline(cin, input);
    timer->TurnOff();
    
    poly->SetPolyLine(-1);
    int npoints = 10;
    for (int i = 0; i < npoints; ++i)
      {
	poly->SetPoint(i, gen.Uniform(0, 35), gen.Uniform(0,35), gen.Uniform(-500, 500));
      }
    poly->Draw();

    c1->Modified();
    c1->Update();
    if (input == "q") done = true;
  } while (!done);
  

}
