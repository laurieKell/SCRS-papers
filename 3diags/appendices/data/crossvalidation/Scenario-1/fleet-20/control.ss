#V3.30.23.2;_safe;_compile_date:_Apr 17 2025;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_13.2
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:_https://groups.google.com/g/ss3-forum_and_NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:_https://nmfs-ost.github.io/ss3-website/
#_Source_code_at:_https://github.com/nmfs-ost/ss3-source-code

#C SS3_Control_NA_SMA_2017_15.xlsx
#C file created using an r4ss function
#C file write time: 2025-06-07  12:54:58
#_data_and_control_files: data.ss // control.ss
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS3)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond sd_ratio_rd < 0: platoon_sd_ratio parameter required after movement params.
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 7 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
1 #_Nblock_Patterns
 1 #_blocks_per_pattern 
# begin and end years of blocks
 1949 1949
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
4 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity;_6=Lorenzen_range
 #_Age_natmort_by sex x growthpattern (nest GP in sex)
 0.206 0.1683 0.1454 0.13 0.119 0.1107 0.1042 0.0991 0.0949 0.0915 0.0887 0.0862 0.0842 0.0824 0.0809 0.0796 0.0785 0.0775 0.0766 0.0758 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751 0.0751
 0.2117 0.1654 0.1423 0.1286 0.1197 0.1136 0.1093 0.1062 0.1039 0.1021 0.1008 0.0998 0.099 0.0984 0.0979 0.0976 0.0973 0.0971 0.0969 0.0967 0.0966 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965 0.0965
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0 #_Age(post-settlement) for L1 (aka Amin); first growth parameter is size at this age; linear growth below this
999 #_Age(post-settlement) for L2 (aka Amax); 999 to treat as Linf
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
4 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
#_Age_Fecundity by growth pattern
 0 0 0 0.01 0.02 0.05 0.12 0.24 0.46 0.79 1.23 1.73 2.26 2.75 3.19 3.55 3.86 4.1 4.3 4.46 4.6 4.71 4.8 4.88 4.95 5.01 5.06 5.1 5.14 5.17 5.2 5.22 5.25 5.27 5.29 5.3 5.32 5.33 5.34 5.35 5.36 5.37 5.37
17 #_First_Mature_Age
# NOTE: maturity options 4 and 5 cause fecundity_at_length to be ignored, but parameters still needed 
2 #_fecundity_at_length option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
# Sex: 1  BioPattern: 1  Growth
 5 100 90.5 90.5 1000 6 -3 0 0 0 0 0.5 0 0 # L_at_Amin_Fem_GP_1
 50 600 361.53 361.53 1000 6 -4 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
 0.01 0.65 0.103 0.103 0.2 6 -5 0 0 0 0 0.5 0 0 # VonBert_K_Fem_GP_1
 0.01 0.3 0.0932677 0.0932677 999 6 -2 0 0 0 0 0.5 0 0 # CV_young_Fem_GP_1
 0.01 0.3 0.15 0.15 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -3 3 5.4e-06 5.4e-06 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Fem_GP_1
 -3 5 3.14 3.14 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 1 300 275.013 275.013 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat50%_Fem_GP_1
 -200 3 -0.10755 -0.10755 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat_slope_Fem_GP_1
 -3 50 6.68e-06 1e-05 0.8 6 -3 0 0 0 0 0.5 0 0 # Eggs_scalar_Fem_GP_1
 -3 3 2.34 2.34442 0.8 6 -3 0 0 0 0 0.5 0 0 # Eggs_exp_len_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
# Sex: 2  BioPattern: 1  Growth
 5 100 87.37 87.37 1000 6 -3 0 0 0 0 0.5 0 0 # L_at_Amin_Mal_GP_1
 50 600 238.49 238.49 1000 6 -4 0 0 0 0 0.5 0 0 # L_at_Amax_Mal_GP_1
 0.01 0.65 0.24 0.24 0.2 6 -5 0 0 0 0 0.5 0 0 # VonBert_K_Mal_GP_1
 0.01 0.3 0.0973418 0.0973418 999 6 -2 0 0 0 0 0.5 0 0 # CV_young_Mal_GP_1
 0.01 0.3 0.15 0.15 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -3 3 1.25e-05 1.25e-05 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Mal_GP_1
 -3 5 2.97 2.97 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
#  Recruitment Distribution 
#  Cohort growth dev base
 0.1 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Platoon StDev Ratio 
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
7 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
           2.3         13.82       5.52809          7.04          1000             6          1          0          0          0          0          0          0          0 # SR_LN(R0)
          0.01             1         0.372       0.39368           0.2             6         -4          0          0          0          0          0          0          0 # SR_surv_zfrac
          0.01            10             3             3           0.2             6         -4          0          0          0          0          0          0          0 # SR_surv_Beta
           0.2           1.9      0.283103          0.28          1000             6         -4          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             6         -4          0          0          0          0          0          0          0 # SR_regime
            -5             5             0             0             1             6         -4          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1990 # first year of main recr_devs; early devs can precede this era
2020 # last year of main recr_devs; forecast devs start in following year
3 #_recdev phase 
1 # (0/1) to read 13 advanced options
 -5 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 4 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1981.59 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1990.48 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2020 #_last_yr_fullbias_adj_in_MPD
 2024 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS3 sets bias_adj to 0.0 for fcast yrs)
 0.3769 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -10 #min rec_dev
 10 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_year Input_value
#
# all recruitment deviations
#  1985E 1986E 1987E 1988E 1989E 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018R 2019R 2020R 2021F 2022F 2023F 2024F
#  -0.150391 -0.159878 -0.166097 -0.162489 -0.134569 -0.180197 -0.108104 0.0698684 -0.0197005 -0.0931695 -0.0329665 -0.236073 -0.355812 -0.348528 -0.0887509 -0.0210717 0.0609049 0.219035 0.0332813 0.291196 0.432976 0.389274 0.00563491 -0.0321536 0.128426 -0.114887 0.0731377 0.0296569 -0.0368514 0.183512 0.0799465 -0.120559 -0.0358415 0.0729557 -0.100718 -0.144422 -0.0441491 0.108083 0.0650494 0
#
#Fishing Mortality info 
0.2 # F ballpark value in units of annual_F
-2010 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
4 # max F (methods 2-4) or harvest fraction (method 1)
5  # N iterations for tuning in hybrid mode; recommend 3 (faster) to 5 (more precise if many fleets)
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 0
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#
# F rates by fleet x season
#_year:  1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# F01_EU_SPN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0241612 0.097963 0.117924 0.0938437 0.0552335 0.078242 0.07067 0.109018 0.0989265 0.108731 0.114481 0.179955 0.138838 0.134784 0.129593 0.0958233 0.0967112 0.109717 0.103854 0.101939 0.0813684 0.0817674 0.0723411 0.0767414 0.0958462 0.0962045 0.0813763 0.11636 0.077461 0.0764867 0.068266 0.0769037 0.088279 0.0581096 0.0419487 0.0412979 0.0322507 0.0317813 0.050614 0.097752
# F02_EU_POR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0108753 0.0176052 0.0123699 0.0446653 0.0360507 0.0374436 0.0417937 0.0224117 0.0204547 0.0229186 0.0218694 0.0241798 0.0246242 0.0693391 0.0251874 0.0569687 0.044907 0.0669766 0.0446956 0.0545307 0.0725761 0.0561258 0.0572771 0.0465739 0.0124074 0.0118941 0.0137857 0.0145201 0.0146416 0.0151494 0.0184384 0.0112267 0.0110195 0.00652971 0.012611
# F03_JPN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00519215 0.00456374 0.00867715 0.00466643 0.00869676 0.00928501 0.00658875 0.0134005 0.0180462 0.00909373 0.0257026 0.0361886 0.0123527 0.0445641 0.00626072 0.00719898 0.00523209 0.0205622 0.0118536 0.0243597 0.0171739 0.0136886 0.00292327 0.00451556 0.00352311 0.00444951 0.00215439 0.00236451 0.00141457 0.00301843 0.00194121 0.00310653 0.00371641 0.00084396 0.00136506 0.0013067 0.000730825 0.000455986 0.000556356 0.0010745
# F04_CTP 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00238564 0.00272616 0.000805403 0.000152027 7.74509e-05 0.000348817 0.00151114 0.000623252 0.00035277 0.00113636 0.00128204 0.00189866 0.00184963 0.00216134 0.00361029 0.00269824 0.00216484 0.00230378 0.00151546 0.00274771 0.00256444 0.00140294 0.000194199 0.000735697 0.000364103 0.000491459 0.000459231 0.000591183 0.000311838 0.000161138 0.000571083 0.000287498 5.17971e-05 0.000952333 0.000208877 0.000489506 5.33397e-05 6.53223e-05 0.000300355 0.000580083
# F05_USA 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00206069 0.0026347 0.00395315 0.00485838 0.00583192 0.00454171 0.00546864 0.00852524 0.00784878 0.00852248 0.00840606 0.00565694 0.00647397 0.00542692 0.00483143 0.00505999 0.0050317 0.00435131 0.00353691 0.00453007 0.00431073 0.00279392 0.00459074 0.0039906 0.00468601 0.00488637 0.00484959 0.00479244 0.0045235 0.00503709 0.00321694 0.00343048 0.00352634 0.00129318 0.00122105 0.00106344 0.00118803 0.00107415 0.00116536 0.00225068
# F06_CAN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00469389 0.00298179 0.00554012 0.00322039 0.00329837 0.00349977 0.00332124 0.00316656 0.0030858 0.00326406 0.00330555 0.0026554 0.00246284 0.00152167 0.00209843 0.00172653 0.00175887 0.001352 0.00172984 0.00265535 0.00409789 0.00391107 0.00522251 0.00299138 0.00310491 0.0019251 0.00184447 0.00227341 0.00075065 0.00144975
# F07_MOR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00818499 0.00898726 0.0110314 0.010398 0.00658417 0.0122391 0.0221928 0.0322601 0.0209538 0.0212853 0.0350971 0.0331027 0.0450098 0.0537201 0.0172232 0.0231062 0.0190805 0.0198102 0.0153906 0.000448752 0.000650249 0.00125584
# F08_MEX 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000111152 0 0 0 0 0.000124278 0.000198953 0 0.000117719 7.89734e-05 0.000112806 6.20821e-05 9.25867e-05 6.74937e-05 7.97028e-05 8.9985e-05 8.21754e-05 8.9588e-05 4.37263e-05 4.17035e-05 4.17973e-05 4.05288e-05 5.65061e-05 2.94022e-05 2.52297e-05 2.74652e-05 2.72872e-05 3.70951e-05 2.45057e-05 4.73284e-05
# F09_VEN 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6.2235e-05 4.0398e-05 6.17442e-05 0.000199946 3.64418e-05 5.33128e-05 1.68515e-05 1.49953e-05 9.10916e-05 0.000112486 0.000325888 9.78833e-05 2.27003e-05 6.07346e-05 6.93273e-05 0.00063822 0.000488592 0.000643726 0.00163425 0.000532433 0.000162249 0.000268002 4.10867e-05 0.000787741 0.000506001 0.000433686 0.000611104 0.000152269 0.00019903 0.000199291 0.000174618 0.00023158 0.000192932 0.000209494 0.000189499 7.07786e-05 1.3338e-05 7.68587e-06 1.48439e-05
# F10_OTH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.131186 0.0605823 0.0226769 0.040319 0.0351317 0.0143375 0.0110978 0.0131893 0.0337454 0.0167958 0.0750802 0.0141362 0.0106906 0.0110965 0.00693246 0.0198709 0.0135726 0.01507 0.00102332 0.0174069 0.014405 0.0116868 0.0125758 0.00975095 0.0133316 0.0117712 0.0162751 0.0224905 0.0225158 0.0464962 0.0326628 0.0121413 0.0211587 0.0197842 0.0117159 0.00226635 0.00183271 0.00189544 0.000539467 0.00104188
# F11_HIS 0.0253663 0.0171255 0.0172198 0.021543 0.00546468 0.0110559 0.00666457 0.0178138 0.0148515 0.0196702 0.0129178 0.0314591 0.0435896 0.0191634 0.0356619 0.0283734 0.0568393 0.0515477 0.0684736 0.0688358 0.0625707 0.100297 0.0959462 0.0922856 0.0992234 0.0764768 0.0831304 0.0898136 0.0942405 0.0933063 0.0743247 0.0982153 0.111351 0.0993768 0.124259 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#
#_Q_setup for fleets with cpue or survey or deviation data
#_1:  fleet number
#_2:  link type: 1=simple q; 2=mirror; 3=power (+1 parm); 4=mirror with scale (+1p); 5=offset (+1p); 6=offset & power (+2p)
#_     where power is applied as y = q * x ^ (1 + power); so a power value of 0 has null effect
#_     and with the offset included it is y = q * (x + offset) ^ (1 + power)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
        12         1         0         0         0         1  #  CPUE_EU_SPN
        13         1         0         0         0         1  #  CPUE_EU_POR
        14         1         0         0         0         1  #  CPUE_JPN_1
        15         1         0         0         0         1  #  CPUE_JPN_2
        16         1         0         0         0         1  #  CPUE_CTP
        17         1         0         0         0         1  #  CPUE_USA
        18         1         0         0         0         1  #  CPUE_MOR_1
        19         1         0         0         0         1  #  CPUE_MOR_2
        20         1         0         0         0         1  #  CPUE_USAlogbook
-9999 0 0 0 0 0
#
#_Q_parameters
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -25            25      -6.42018             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_EU_SPN(12)
           -25            25      -9.83094             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_EU_POR(13)
           -25            25      -6.21484             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_JPN_1(14)
           -25            25      -6.29811             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_JPN_2(15)
           -25            25      -6.46477             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_CTP(16)
           -25            25       -6.5631             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_USA(17)
           -25            25      -9.65124             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_MOR_1(18)
           -25            25      -11.6164             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_MOR_2(19)
           -25            25      -6.92571             0          0.05             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_USAlogbook(20)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0;  parm=0; selex=1.0 for all sizes
#Pattern:_1;  parm=2; logistic; with 95% width specification
#Pattern:_5;  parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_11; parm=2; selex=1.0  for specified min-max population length bin range
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6;  parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (mean over bin range)
#Pattern:_8;  parm=8; double_logistic with smooth transitions and constant above Linf option
#Pattern:_9;  parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2*special; non-parm len selex, read as N break points, then N selex parameters
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_2;  parm=6; double_normal with sel(minL) and sel(maxL), using joiners, back compatibile version of 24 with 3.30.18 and older
#Pattern:_25; parm=3; exponential-logistic in length
#Pattern:_27; parm=special+3; cubic spline in length; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=special+3+2; cubic spline; like 27, with 2 additional param for scaling (mean over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 0 0 0 # 1 F01_EU_SPN
 24 0 0 0 # 2 F02_EU_POR
 24 0 3 0 # 3 F03_JPN
 24 0 3 0 # 4 F04_CTP
 24 0 4 0 # 5 F05_USA
 24 0 0 0 # 6 F06_CAN
 15 0 0 2 # 7 F07_MOR
 24 0 3 0 # 8 F08_MEX
 24 0 0 0 # 9 F09_VEN
 15 0 0 1 # 10 F10_OTH
 15 0 0 1 # 11 F11_HIS
 15 0 0 1 # 12 CPUE_EU_SPN
 15 0 0 2 # 13 CPUE_EU_POR
 15 0 0 3 # 14 CPUE_JPN_1
 15 0 0 3 # 15 CPUE_JPN_2
 15 0 0 4 # 16 CPUE_CTP
 15 0 0 5 # 17 CPUE_USA
 15 0 0 7 # 18 CPUE_MOR_1
 15 0 0 7 # 19 CPUE_MOR_2
 15 0 0 5 # 20 CPUE_USAlogbook
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic. Recommend using pattern 18 instead.
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (mean over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (mean over bin range)
#Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
 11 0 0 0 # 1 F01_EU_SPN
 15 0 0 1 # 2 F02_EU_POR
 15 0 0 1 # 3 F03_JPN
 15 0 0 1 # 4 F04_CTP
 15 0 0 1 # 5 F05_USA
 15 0 0 1 # 6 F06_CAN
 15 0 0 1 # 7 F07_MOR
 15 0 0 1 # 8 F08_MEX
 15 0 0 1 # 9 F09_VEN
 15 0 0 1 # 10 F10_OTH
 15 0 0 1 # 11 F11_HIS
 15 0 0 1 # 12 CPUE_EU_SPN
 15 0 0 1 # 13 CPUE_EU_POR
 15 0 0 1 # 14 CPUE_JPN_1
 15 0 0 1 # 15 CPUE_JPN_2
 15 0 0 1 # 16 CPUE_CTP
 15 0 0 1 # 17 CPUE_USA
 15 0 0 1 # 18 CPUE_MOR_1
 15 0 0 1 # 19 CPUE_MOR_2
 15 0 0 1 # 20 CPUE_USAlogbook
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   F01_EU_SPN LenSelex
          62.5         297.5       123.757        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F01_EU_SPN(1)
           -10             4      -4.36683            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F01_EU_SPN(1)
            -1             9       6.10381           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F01_EU_SPN(1)
            -1             9       7.87496          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F01_EU_SPN(1)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F01_EU_SPN(1)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F01_EU_SPN(1)
# 2   F02_EU_POR LenSelex
          62.5         297.5       141.641        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F02_EU_POR(2)
           -10             4      -7.80225            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F02_EU_POR(2)
            -1             9       6.68926           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F02_EU_POR(2)
            -1             9       7.17732          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F02_EU_POR(2)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F02_EU_POR(2)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F02_EU_POR(2)
# 3   F03_JPN LenSelex
          62.5         297.5       150.103        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F03_JPN(3)
           -10             4      -7.28993            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F03_JPN(3)
            -1             9       7.17117           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F03_JPN(3)
            -1             9       7.32804          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F03_JPN(3)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F03_JPN(3)
            -5             9         -4.96            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F03_JPN(3)
           -60            60      -1.80654             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Peak_F03_JPN(3)
           -15            15     0.0633349             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Ascend_F03_JPN(3)
           -15            15       0.50468             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Descend_F03_JPN(3)
           -15            15             0             0          0.05             1         -4          0          0          0          0          0          0          0  #  SzSel_Male_Final_F03_JPN(3)
           0.5           1.5      0.828244             1          0.05             1          5          0          0          0          0          0          0          0  #  SzSel_Male_Scale_F03_JPN(3)
# 4   F04_CTP LenSelex
          62.5         297.5        127.28        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F04_CTP(4)
            -6             4      -2.19096            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F04_CTP(4)
            -1             9       3.75874           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F04_CTP(4)
            -1             9       7.83374          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F04_CTP(4)
           -10             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F04_CTP(4)
           -10             9         -4.96            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F04_CTP(4)
           -60            60     0.0718952             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Peak_F04_CTP(4)
           -15            15       1.42188             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Ascend_F04_CTP(4)
           -15            15     -0.491315             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Descend_F04_CTP(4)
           -15            15             0             0          0.05             1         -4          0          0          0          0          0          0          0  #  SzSel_Male_Final_F04_CTP(4)
           0.5           1.5      0.905123             1          0.05             1          5          0          0          0          0          0          0          0  #  SzSel_Male_Scale_F04_CTP(4)
# 5   F05_USA LenSelex
          62.5         297.5       134.048        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F05_USA(5)
           -10             4      -2.36144            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F05_USA(5)
            -1            15       4.85897           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F05_USA(5)
            -1             9       7.76278          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F05_USA(5)
           -10             9     -0.519001            -5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F05_USA(5)
           -10             9      -2.86875            -5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F05_USA(5)
           -40            40        -10.95             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F05_USA(5)
           -15            15      -7.61495             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F05_USA(5)
           -15            15       0.11646             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F05_USA(5)
           -15            15      -2.38976             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F05_USA(5)
           0.5           1.5      0.957881             1          0.05             1          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F05_USA(5)
# 6   F06_CAN LenSelex
          62.5         297.5       133.933        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F06_CAN(6)
           -10             4      -6.54961            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F06_CAN(6)
            -1             9       7.21324           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F06_CAN(6)
            -1             9       7.62939          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F06_CAN(6)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F06_CAN(6)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F06_CAN(6)
# 7   F07_MOR LenSelex
# 8   F08_MEX LenSelex
          62.5         297.5       189.945           188          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F08_MEX(8)
            -6            10      -5.27113          -5.3          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F08_MEX(8)
            -1            10       8.21579             8          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F08_MEX(8)
            -1            10       8.56406             8          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F08_MEX(8)
           -10             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F08_MEX(8)
           -10            10      -2.61098            -5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F08_MEX(8)
           -60            60       -20.255             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Peak_F08_MEX(8)
           -15            15      -0.36768             0          0.05             1          4          0          0          0          0          0          0          0  #  SzSel_Male_Ascend_F08_MEX(8)
           -15            15            15            15          0.05             1         -4          0          0          0          0          0          0          0  #  SzSel_Male_Descend_F08_MEX(8)
           -15            15            15            15          0.05             1         -4          0          0          0          0          0          0          0  #  SzSel_Male_Final_F08_MEX(8)
           0.5           1.5      0.909146             1          0.05             1          5          0          0          0          0          0          0          0  #  SzSel_Male_Scale_F08_MEX(8)
# 9   F09_VEN LenSelex
          62.5         297.5       177.397        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F09_VEN(9)
            -6             4      -5.30633            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F09_VEN(9)
            -1             9       7.65844           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F09_VEN(9)
            -1             9       7.43349          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F09_VEN(9)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F09_VEN(9)
           -10             9         -4.96            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F09_VEN(9)
# 10   F10_OTH LenSelex
# 11   F11_HIS LenSelex
# 12   CPUE_EU_SPN LenSelex
# 13   CPUE_EU_POR LenSelex
# 14   CPUE_JPN_1 LenSelex
# 15   CPUE_JPN_2 LenSelex
# 16   CPUE_CTP LenSelex
# 17   CPUE_USA LenSelex
# 18   CPUE_MOR_1 LenSelex
# 19   CPUE_MOR_2 LenSelex
# 20   CPUE_USAlogbook LenSelex
# 1   F01_EU_SPN AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F01_EU_SPN(1)
            10           100            42             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F01_EU_SPN(1)
# 2   F02_EU_POR AgeSelex
# 3   F03_JPN AgeSelex
# 4   F04_CTP AgeSelex
# 5   F05_USA AgeSelex
# 6   F06_CAN AgeSelex
# 7   F07_MOR AgeSelex
# 8   F08_MEX AgeSelex
# 9   F09_VEN AgeSelex
# 10   F10_OTH AgeSelex
# 11   F11_HIS AgeSelex
# 12   CPUE_EU_SPN AgeSelex
# 13   CPUE_EU_POR AgeSelex
# 14   CPUE_JPN_1 AgeSelex
# 15   CPUE_JPN_2 AgeSelex
# 16   CPUE_CTP AgeSelex
# 17   CPUE_USA AgeSelex
# 18   CPUE_MOR_1 AgeSelex
# 19   CPUE_MOR_2 AgeSelex
# 20   CPUE_USAlogbook AgeSelex
#_No_Dirichlet parameters
#_no timevary selex parameters
#
0   #  use 2D_AR1 selectivity? (0/1)
#_no 2D_AR1 selex offset used
#_specs:  fleet, ymin, ymax, amin, amax, sigma_amax, use_rho, len1/age2, devphase, before_range, after_range
#_sigma_amax>amin means create sigma parm for each bin from min to sigma_amax; sigma_amax<0 means just one sigma parm is read and used for all bins
#_needed parameters follow each fleet's specifications
# -9999  0 0 0 0 0 0 0 0 0 0 # terminator
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# no timevary parameters
#
#
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_factor  fleet  value
      1     12     0.132
      1     13         0
      1     14     0.163
      1     15         0
      1     16     0.165
      1     17     0.022
      1     18     0.129
      1     19      0.01
      1     20      0.12
      4      1  0.455705
      4      2  0.639037
      4      3   1.01233
      4      4  0.321861
      4      5   1.51484
      4      6  0.441869
      4      7   1.00625
      4      8   3.07597
      4      9   4.60254
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 9 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 1 12 1 1 1
 1 13 1 1 1
 1 14 1 1 1
 1 15 1 1 1
 1 16 1 1 1
 1 17 1 1 1
 1 18 1 1 1
 1 19 1 1 1
 1 20 1 0 1
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  0 #_CPUE/survey:_3
#  0 #_CPUE/survey:_4
#  0 #_CPUE/survey:_5
#  0 #_CPUE/survey:_6
#  0 #_CPUE/survey:_7
#  0 #_CPUE/survey:_8
#  0 #_CPUE/survey:_9
#  0 #_CPUE/survey:_10
#  0 #_CPUE/survey:_11
#  1 #_CPUE/survey:_12
#  1 #_CPUE/survey:_13
#  1 #_CPUE/survey:_14
#  1 #_CPUE/survey:_15
#  1 #_CPUE/survey:_16
#  1 #_CPUE/survey:_17
#  1 #_CPUE/survey:_18
#  1 #_CPUE/survey:_19
#  0 #_CPUE/survey:_20
#  1 #_lencomp:_1
#  1 #_lencomp:_2
#  1 #_lencomp:_3
#  1 #_lencomp:_4
#  1 #_lencomp:_5
#  1 #_lencomp:_6
#  1 #_lencomp:_7
#  1 #_lencomp:_8
#  1 #_lencomp:_9
#  0 #_lencomp:_10
#  0 #_lencomp:_11
#  0 #_lencomp:_12
#  0 #_lencomp:_13
#  0 #_lencomp:_14
#  0 #_lencomp:_15
#  0 #_lencomp:_16
#  0 #_lencomp:_17
#  0 #_lencomp:_18
#  0 #_lencomp:_19
#  0 #_lencomp:_20
#  1 #_init_equ_catch1
#  1 #_init_equ_catch2
#  1 #_init_equ_catch3
#  1 #_init_equ_catch4
#  1 #_init_equ_catch5
#  1 #_init_equ_catch6
#  1 #_init_equ_catch7
#  1 #_init_equ_catch8
#  1 #_init_equ_catch9
#  1 #_init_equ_catch10
#  1 #_init_equ_catch11
#  1 #_init_equ_catch12
#  1 #_init_equ_catch13
#  1 #_init_equ_catch14
#  1 #_init_equ_catch15
#  1 #_init_equ_catch16
#  1 #_init_equ_catch17
#  1 #_init_equ_catch18
#  1 #_init_equ_catch19
#  1 #_init_equ_catch20
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  0 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999

