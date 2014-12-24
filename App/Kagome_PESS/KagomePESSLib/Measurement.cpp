
double measureObs(map<string, UniTensor> &La, map<string, UniTensor> &groupT , string Tname, UniTensor &Obs, string Obsname, string PATH);
double measureNorm(UniTensor &T);
double Measurement(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, UniTensor &Obs);

double measureObs(map<string, UniTensor> &La, map<string, UniTensor> &US , string Tname, UniTensor &Obs, string Obsname, string PATH){

  for(map<string, UniTensor>::iterator ita = La.begin(); ita != La.end(); ++ita){
    string U = ita -> first;
    U[0] = 'U';
    bondcat(US[U], ita->second.getBlock(), 1);
  }

//  UniTensor T = TLaunch(groupT);
  Network contractTup("M.net");
  contractTup.putTensor("U0", &US["U0"]);
  contractTup.putTensor("U1", &US["U1"]);
  contractTup.putTensor("U2", &US["U2"]);
  contractTup.putTensor("S0", &US["S0"]);

  UniTensor T = contractTup.launch();
  Network measureobs("KagomePESSnet/MeasureH_T.net");
  measureobs.putTensor("T", T);
  measureobs.putTensor("TT", T);
  measureobs.putTensor("H", Obs);


  /*
  map<string, UniTensor> Measure;
  Measure[Tname] = T;
  Measure[Obsname] = Obs;
  char TTname[8];
  strcpy(TTname, Tname.c_str());
  strcat(TTname, "T");
  Measure[TTname] = T;

  //Open .netfile
  char netfile[32];
  char buf[32];
  strcpy(netfile,"Measure");
  sprintf(buf, "%s_%s.net", Obsname.c_str(), Tname.c_str());
  strcat(netfile, buf);
  string filedir;
  filedir.assign(netfile);

  FILE *fp;
  if((fp = fopen(netFileDir(PATH, filedir).c_str(), "r")) == NULL)
    writeMeasureNetFile(Measure, Tname, Obsname, netFileDir(PATH, filedir));
  else
    fclose(fp);

  Network measureobs(netFileDir(PATH, filedir));
  for(map<string, UniTensor>::iterator it = Measure.begin(); it != Measure.end(); ++it)
    measureobs.putTensor(it->first, it->second);

    */
  double val = measureobs.launch()[0];

  double MN = (T*T)[0];//measureNorm(T);
  double measure = val/MN;

  return measure;
}

double measureNorm(UniTensor &T){

  Network norm("Norm.net");
  norm.putTensor("T", &T);
  norm.putTensorT("TT", &T);
  UniTensor M = norm.launch();

  double measurenorm = M[0];

  return measurenorm;
}

double Measurement(map<string, UniTensor> &groupTu, map<string, UniTensor> &groupTd, UniTensor &Obs){

  Network MeasurementT("MeasurementT.net");
  MeasurementT.putTensor("U0", &groupTu["U0"]);
  MeasurementT.putTensor("U1", &groupTu["U1"]);
  MeasurementT.putTensor("U2", &groupTu["U2"]);
  MeasurementT.putTensor("S0u", &groupTu["S0"]);
  MeasurementT.putTensor("U3", &groupTd["U1"]);
  MeasurementT.putTensor("U4", &groupTd["U2"]);
  MeasurementT.putTensor("S0d", &groupTd["S0"]);
  UniTensor T = MeasurementT.launch();

  Network MeasurementObs("MeasurementObs.net");
  MeasurementObs.putTensor("T", &T);
  MeasurementObs.putTensor("TT", &T);
  MeasurementObs.putTensor("Obs", &Obs);
  UniTensor M = MeasurementObs.launch();

  Network MeasurementNorm("MeasurementNorm.net");
  MeasurementNorm.putTensor("T", &T);
  MeasurementNorm.putTensor("TT", &T);
  UniTensor N = MeasurementNorm.launch();

  double Measure = M[0] / N[0];
  //  assert(groupTu["U2"] == groupTd["U0"]);
  return Measure;
}
