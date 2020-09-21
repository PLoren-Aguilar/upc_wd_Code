      MODULE mod_EOS
!===============================================
! declaration for pipelining the eos routines
!===============================================

! maximum length of the row vector
      integer, parameter :: nrowmax = 1000

! maximum number of isotopes
      integer, parameter :: irowmax = 20

! maximum number of ionization stages
      integer, parameter :: jstagemax = 1

! io variables for helmeos
      integer :: jlo_eos, jhi_eos
      double precision, dimension(nrowmax) :: ewant_row, temp_row,      &
      den_row, abar_row, zbar_row, det_row, ptot_row, cs_row, cv_row,   &
      dpt_row, cp_row, etot_row
      integer, dimension(nrowmax) :: eoslist
      logical :: eosfail

! thermodynamic and composition inputs
      double precision, dimension(nrowmax) :: zeff_row, ye_row

! composition input
      integer :: niso
      double precision, dimension(irowmax,nrowmax) :: xmass_row

! composition output
      double precision, dimension(jstagemax,irowmax,nrowmax) :: frac_row

! composition output for sneos
      double precision, dimension(nrowmax) :: xn_row, xp_row, xa_row,   &
      xhv_row, xmuhat_row

! totals and their derivatives
      double precision, dimension(nrowmax) :: dpd_row, dpa_row,dpz_row, &
      dpdd_row, dpdt_row, dpda_row, dpdz_row, dptt_row, dpta_row,       &
      dptz_row, dpaa_row, dpaz_row, dpzz_row

      double precision, dimension(nrowmax) :: ded_row, dea_row, dez_row,&
      dedd_row, dedt_row, deda_row, dedz_row, dett_row, deta_row,       &
      detz_row, deaa_row, deaz_row, dezz_row

      double precision, dimension(nrowmax) :: stot_row, dst_row,        &
      dsd_row, dsa_row, dsz_row, dsdd_row, dsdt_row, dsda_row, dsdz_row,&
      dstt_row, dsta_row, dstz_row, dsaa_row, dsaz_row, dszz_row

! radiation contributions
      double precision, dimension(nrowmax) :: prad_row, dpradt_row,     &
      dpradd_row, dprada_row, dpradz_row, dpraddd_row, dpraddt_row,     &
      dpradda_row, dpraddz_row, dpradtt_row, dpradta_row, dpradtz_row,  &
      dpradaa_row, dpradaz_row, dpradzz_row

      double precision, dimension(nrowmax) :: erad_row, deradt_row,     &
      deradd_row, derada_row, deradz_row, deraddd_row, deraddt_row,     &
      deradda_row, deraddz_row, deradtt_row, deradta_row, deradtz_row,  &
      deradaa_row, deradaz_row, deradzz_row

      double precision, dimension(nrowmax) :: srad_row, dsradt_row,     &
      dsradd_row, dsrada_row, dsradz_row, dsraddd_row, dsraddt_row,     &
      dsradda_row, dsraddz_row, dsradtt_row, dsradta_row, dsradtz_row,  &
      dsradaa_row, dsradaz_row, dsradzz_row

! gas contributions
      double precision, dimension(nrowmax) :: pgas_row, dpgast_row,     &
      dpgasd_row, dpgasa_row, dpgasz_row, dpgasdd_row, dpgasdt_row,     &
      dpgasda_row, dpgasdz_row, dpgastt_row, dpgasta_row, dpgastz_row,  &
      dpgasaa_row, dpgasaz_row, dpgaszz_row

      double precision, dimension(nrowmax) :: egas_row, degast_row,     &
      degasd_row, degasa_row, degasz_row, degasdd_row, degasdt_row,     &
      degasda_row, degasdz_row, degastt_row, degasta_row, degastz_row,  &
      degasaa_row, degasaz_row, degaszz_row

      double precision, dimension(nrowmax) :: sgas_row, dsgast_row,     &
      dsgasd_row, dsgasa_row, dsgasz_row, dsgasdd_row, dsgasdt_row,     &
      dsgasda_row, dsgasdz_row, dsgastt_row, dsgasta_row, dsgastz_row,  &
      dsgasaa_row, dsgasaz_row, dsgaszz_row

! ion contributions
      double precision, dimension(nrowmax) :: pion_row, dpiont_row,     &
      dpiond_row, dpiona_row, dpionz_row, dpiondd_row, dpiondt_row,     &
      dpionda_row, dpiondz_row, dpiontt_row, dpionta_row, dpiontz_row,  &
      dpionaa_row, dpionaz_row, dpionzz_row

      double precision, dimension(nrowmax) :: eion_row, deiont_row,     &
      deiond_row, deiona_row,deionz_row, deiondd_row, deiondt_row,      &
      deionda_row, deiondz_row, deiontt_row, deionta_row, deiontz_row,  &
      deionaa_row, deionaz_row, deionzz_row

      double precision, dimension(nrowmax) :: sion_row, dsiont_row,     &
      dsiond_row, dsiona_row, dsionz_row, dsiondd_row, dsiondt_row,     &
      dsionda_row, dsiondz_row, dsiontt_row, dsionta_row, dsiontz_row,  &
      dsionaa_row, dsionaz_row, dsionzz_row

      double precision, dimension(nrowmax) :: etaion_row, detait_row,   &
      detaid_row, detaia_row, detaiz_row, detaidd_row, detaidt_row,     &
      detaida_row, detaidz_row, detaitt_row, detaita_row, detaitz_row,  &
      detaiaa_row, detaiaz_row, detaizz_row

      double precision, dimension(nrowmax) :: xni_row, xnim_row,        &
      dxnit_row, dxnid_row, dxnia_row, dxniz_row, dxnidd_row,dxnidt_row,&
      dxnida_row, dxnidz_row, dxnitt_row, dxnita_row, dxnitz_row,       &
      dxniaa_row, dxniaz_row, dxnizz_row

! electron-positron contributions
      double precision, dimension(nrowmax) :: etaele_row, etapos_row,   &
      detat_row, detad_row, detaa_row, detaz_row, detadd_row,           &
      detadt_row, detada_row, detadz_row, detatt_row, detata_row,       &
      detatz_row, detaaa_row, detaaz_row, detazz_row

      double precision, dimension(nrowmax) :: pele_row, ppos_row,       &
      dpept_row, dpepd_row, dpepa_row, dpepz_row, dpepdd_row,dpepdt_row,&
      dpepda_row, dpepdz_row, dpeptt_row, dpepta_row, dpeptz_row,       &
      dpepaa_row, dpepaz_row, dpepzz_row

      double precision, dimension(nrowmax) :: eele_row, epos_row,       &
      deept_row, deepd_row, deepa_row, deepz_row, deepdd_row,deepdt_row,&
      deepda_row, deepdz_row, deeptt_row, deepta_row, deeptz_row,       &
      deepaa_row, deepaz_row, deepzz_row

      double precision, dimension(nrowmax) :: sele_row, spos_row,        &
      dsept_row, dsepd_row, dsepa_row, dsepz_row, dsepdd_row, dsepdt_row,&
      dsepda_row, dsepdz_row, dseptt_row, dsepta_row, dseptz_row,        &
      dsepaa_row, dsepaz_row, dsepzz_row

      double precision, dimension(nrowmax) :: xne_row, xnp_row, xnem_row, &
      dxnet_row, dxned_row, dxnea_row, dxnez_row, dxnedd_row, dxnedt_row, &
      dxneda_row, dxnedz_row, dxnett_row, dxneta_row, dxnetz_row,         &
      dxneaa_row, dxneaz_row, dxnezz_row

! ionization potential contributions
      double precision, dimension(nrowmax) :: pip_row, eip_row, sip_row

! coulomb contributions
      double precision, dimension(nrowmax) :: pcou_row, dpcout_row,     &
      dpcoud_row, dpcoua_row, dpcouz_row, ecou_row, decout_row,         &
      decoud_row, decoua_row, decouz_row, scou_row, dscout_row,         &
      dscoud_row, dscoua_row, dscouz_row, plasg_row

! thermodynamic consistency checks; maxwell relations
      double precision, dimension(nrowmax) :: dse_row, dpe_row, dsp_row

! derivative based quantities for the gas
      double precision, dimension(nrowmax) :: cp_gas_row, dcp_gasdd_row,&
      dcp_gasdt_row, dcp_gasda_row, dcp_gasdz_row, cv_gas_row,          &
      dcv_gasdd_row, dcv_gasdt_row, dcv_gasda_row, dcv_gasdz_row

      double precision, dimension(nrowmax) :: gam1_gas_row, dgam1_gasdd_row, &
      dgam1_gasdt_row, dgam1_gasda_row, dgam1_gasdz_row, gam2_gas_row,       &
      dgam2_gasdd_row, dgam2_gasdt_row, dgam2_gasda_row, dgam2_gasdz_row,    &
      gam3_gas_row, dgam3_gasdd_row, dgam3_gasdt_row, dgam3_gasda_row,       &
      dgam3_gasdz_row, nabad_gas_row, dnab_gasdd_row, dnab_gasdt_row,        &
      dnab_gasda_row, dnab_gasdz_row, cs_gas_row, dcs_gasdd_row,             &
      dcs_gasdt_row, dcs_gasda_row, dcs_gasdz_row

! derivative based quantities for the totals
      double precision, dimension(nrowmax) :: dcpdd_row, dcpdt_row,     &
      dcpda_row, dcpdz_row, dcvdd_row, dcvdt_row, dcvda_row, dcvdz_row

      double precision, dimension(nrowmax) :: gam1_row, dgam1dd_row,    &
      dgam1dt_row, dgam1da_row, dgam1dz_row, gam2_row, dgam2dd_row,     &
      dgam2dt_row, dgam2da_row, dgam2dz_row, gam3_row, dgam3dd_row,     &
      dgam3dt_row, dgam3da_row, dgam3dz_row, nabad_row, dnabdd_row,     &
      dnabdt_row, dnabda_row, dnabdz_row, dcsdd_row, dcsdt_row,         &
      dcsda_row, dcsdz_row

! a few work arrays
      double precision, dimension(nrowmax) :: eoswrk01, eoswrk02,       &
      eoswrk03, eoswrk04

! for debugging
      double precision, dimension(nrowmax) :: crp_row, dcrpt_row,         &
      dcrpd_row, dcrpa_row,dcrpz_row, dcrpdd_row, dcrpdt_row, dcrpda_row, &
      dcrpdz_row, dcrptt_row, dcrpta_row, dcrptz_row, dcrpaa_row,         &
      dcrpaz_row, dcrpzz_row
!
      END MODULE mod_EOS
