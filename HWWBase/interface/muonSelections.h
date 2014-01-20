#ifndef MUON_SELECTIONS_H
#define MUON_SELECTIONS_H

namespace HWWFunctions {

  ///////////////
  // Selectors //
  ///////////////

  enum SelectionType { 

      ///////////////
      // Higgs, WW //
      ///////////////

        // WW
        NominalWWV0,
        muonSelectionFO_mu_ww,
        muonSelectionFO_mu_ww_iso10,
        NominalWWV1,
        muonSelectionFO_mu_wwV1,
        muonSelectionFO_mu_wwV1_iso10,
        muonSelectionFO_mu_wwV1_iso10_d0,

        // Analysis
        NominalSmurfV3,
        NominalSmurfV4,
        NominalSmurfV5,
        NominalSmurfV6,

        // Fakes
        muonSelectionFO_mu_smurf_04,
        muonSelectionFO_mu_smurf_10,
  }; 

  ////////////////////
  // Identification //
  ////////////////////

  bool muonId           ( unsigned int index, SelectionType type);
  bool muonIdNotIsolated( unsigned int index, SelectionType type);
  bool isGoodStandardMuon( unsigned int index );

  ///////////////
  // Isolation //
  ///////////////
  double muonIsoValue          ( unsigned int , bool = true );
  double muonIsoValue_TRK      ( unsigned int , bool = true );
  double muonIsoValue_ECAL     ( unsigned int , bool = true );
  double muonIsoValue_HCAL     ( unsigned int , bool = true );

  double muonIsoValuePF        ( unsigned int imu, unsigned int ivtx, float coner=0.4, float minptn=1.0, float dzcut=0.1, int filterId = 0);

  /////////////////////////////
  // Muon d0 corrected by PV //
  /////////////////////////////

  double mud0PV         (unsigned int index);
  double mud0PV_wwV1    (unsigned int index);
  double mudzPV_wwV1    (unsigned int index);
  double mud0PV_smurfV3 (unsigned int index);
  double mudzPV_smurfV3 (unsigned int index);

  /////////////////////////////////////////////////////////////////
  // checks if muon is also pfmuon, and pfmuon pt = reco muon pt //
  /////////////////////////////////////////////////////////////////

  bool isPFMuon( int index , bool requireSamePt = true , float dpt_max = 1.0 );

  struct mu2012_tightness 
  {
      enum value_type 
      {
          LOOSE,
          TIGHT,
          static_size
      };
  };

}
#endif
