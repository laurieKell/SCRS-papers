#V3.30.23.2;_safe;_compile_date:_Apr 17 2025;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_13.2
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:_https://groups.google.com/g/ss3-forum_and_NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:_https://nmfs-ost.github.io/ss3-website/
#_Source_code_at:_https://github.com/nmfs-ost/ss3-source-code

#C SS3_Control_NA_SMA_2017_15.xlsx																																																																																																																																																																																								
#_data_and_control_files: data.ss // control.ss
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS3)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond sd_ratio_rd < 0: platoon_sd_ratio parameter required after movement params.
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
2 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
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
 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0795 0.0785 0.0777 0.0769 0.0762 0.0749 #natm
 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0801 0.0795 0.0785 0.0777 0.0769 0.0762 0.0749 #natm
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
 0 0 0 0 0 0 0 0 0 0 0 0 0.01 0.01 0.03 0.08 0.16 0.3 0.54 0.88 1.32 1.81 2.29 2.74 3.13 3.45 3.71 3.93 4.12 4.28 4.42
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
 5 100 62.9093 62.9093 1000 6 -3 0 0 0 0 0.5 0 0 # L_at_Amin_Fem_GP_1
 50 600 350.6 350.6 1000 6 -4 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
 0.01 0.65 0.064 0.064 0.2 6 -5 0 0 0 0 0.5 0 0 # VonBert_K_Fem_GP_1
 0.01 0.3 0.0932677 0.0932677 999 6 -2 0 0 0 0 0.5 0 0 # CV_young_Fem_GP_1
 0.01 0.3 0.0898002 0.0898002 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -3 3 5.24e-06 5.24e-06 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Fem_GP_1
 -3 5 3.14 3.14 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 1 300 275.013 275.013 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat50%_Fem_GP_1
 -200 3 -0.10755 -0.10755 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat_slope_Fem_GP_1
 -3 50 6.68e-06 1e-05 0.8 6 -3 0 0 0 0 0.5 0 0 # Eggs_scalar_Fem_GP_1
 -3 3 2.34 2.34442 0.8 6 -3 0 0 0 0 0.5 0 0 # Eggs_exp_len_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
# Sex: 2  BioPattern: 1  Growth
 5 100 63.0138 63.0138 1000 6 -3 0 0 0 0 0.5 0 0 # L_at_Amin_Mal_GP_1
 50 600 241.8 241.8 1000 6 -4 0 0 0 0 0.5 0 0 # L_at_Amax_Mal_GP_1
 0.01 0.65 0.136 0.136 0.2 6 -5 0 0 0 0 0.5 0 0 # VonBert_K_Mal_GP_1
 0.01 0.3 0.0973418 0.0973418 999 6 -2 0 0 0 0 0.5 0 0 # CV_young_Mal_GP_1
 0.01 0.3 0.0824247 0.0824247 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -3 3 5.24e-06 5.24e-06 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Mal_GP_1
 -3 5 3.14 3.14 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
#  Recruitment Distribution 
 -4 4 0 1 99 0 -3 0 0 0 0 0.5 0 0 # RecrDist_GP_1
 -4 4 0 1 99 0 -3 0 0 0 0 0.5 0 0 # RecrDist_Area_1
 -4 4 0 1 99 0 -3 0 0 0 0 0.5 0 0 # RecrDist_month_1
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
           2.3         13.82       5.39769          7.04          1000             6          1          0          0          0          0          0          0          0 # SR_LN(R0)
          0.01             1         0.171           0.5           0.2             6         -4          0          0          0          0          0          0          0 # SR_surv_zfrac
          0.01            10             3             1           0.2             6         -4          0          0          0          0          0          0          0 # SR_surv_Beta
           0.2           1.9      0.283103          0.28          1000             6         -4          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             6         -4          0          0          0          0          0          0          0 # SR_regime
            -5             5             0             0             1             6         -4          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1990 # first year of main recr_devs; early devs can precede this era
2012 # last year of main recr_devs; forecast devs start in following year
3 #_recdev phase 
1 # (0/1) to read 13 advanced options
 -5 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 4 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1981.59 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1990.48 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2012 #_last_yr_fullbias_adj_in_MPD
 2016 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS3 sets bias_adj to 0.0 for fcast yrs)
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
#  1985E 1986E 1987E 1988E 1989E 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013F 2014F 2015F 2016F
#  -0.180378 -0.247134 -0.260448 -0.242128 -0.149023 -0.181066 -0.196887 -0.247507 -0.251714 -0.288855 -0.392899 -0.394369 -0.172802 0.0870386 0.0675307 0.253124 0.256116 -0.0539 0.348977 0.460551 0.483154 0.213748 -0.0344924 0.0652315 0.193521 -0.0299442 -0.100972 -0.0835844 -0.0374922 -0.031738 -0.0436861 0
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
#_year:  1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# F1_EU_LL 0.00358033 0.00239757 0.00240041 0.00299261 0.000759044 0.00154006 0.000930254 0.0024872 0.00207009 0.00273809 0.00179781 0.00423606 0.00574182 0.00249959 0.00450199 0.00358693 0.00751724 0.00675742 0.0089409 0.0088407 0.00799285 0.00857776 0.00815654 0.00975779 0.00738321 0.00959793 0.00722484 0.00849074 0.00927407 0.006634 0.00980941 0.0104055 0.0119372 0.02178 0.0142845 0.0204541 0.0826939 0.100468 0.0823327 0.0506471 0.0851195 0.0873912 0.126338 0.15602 0.16942 0.186717 0.289646 0.222502 0.21745 0.20963 0.158382 0.159497 0.174423 0.220584 0.159541 0.176889 0.167768 0.185 0.155665 0.182424 0.204025 0.169532 0.224785 0.171549 0.13348 0.132406 0.0266384
# F2_JPN_LL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00250832 0.00258342 0.00137387 0.00693649 0.00780691 0.00191163 0.00537742 0.00349218 0.00102805 0.00563165 0.00891972 0.00636179 0.00374931 0.00336843 0.00348278 0.00307125 0.00585407 0.00319005 0.00612187 0.00682883 0.005076 0.010795 0.0153618 0.00828208 0.0248329 0.036825 0.0132887 0.0498147 0.00711677 0.0083337 0.00622327 0.0250351 0.0147032 0.0304962 0 0 0.00387026 0.00587931 0.00432967 0.00525122 0.0025124 0.00279094 0.001742 0.00391065 0.00281765 0.000566875
# F3_CTP_LL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000703407 0.00115468 0.00132506 0.00159242 0.00165549 0.001892 0.000558171 0.000106385 5.56049e-05 0.000261083 0.00118375 0.000510668 0.000306411 0.00106051 0.00127122 0.00198514 0.00205159 0.00250454 0.00429393 0.00333075 0.00281381 0.00313295 0.0021392 0.00394787 0.00370822 0.00208328 0.000297262 0.00127399 0.00406627 0.00064333 0.00253681 0.00170286 0.000667025 0.00086121 0.000644434 0.000129652
# F4_USA_LL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000776597 0.000787467 0.000802469 0.00100636 0.00129269 0.00182211 0.0023506 0.00284881 0.00224682 0.00284028 0.00424454 0.00846388 0.00994196 0.010037 0.00839363 0.00954005 0.0082872 0.00402291 0.00750525 0.00817511 0.00738484 0.00611149 0.00792006 0.0076476 0.00509091 0.00841819 0.00718344 0.00795836 0.00833374 0.00816736 0.00798249 0.00812644 0.00931569 0.0161775 0.0032547
# F5_VEN_LL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.81489e-05 3.13626e-05 4.82131e-05 0.000158408 2.9648e-05 4.4906e-05 1.47414e-05 1.36666e-05 8.74467e-05 0.000114372 0.000348728 0.000110513 2.69613e-05 7.48346e-05 8.81647e-05 0.000843278 0.000674761 0.000926604 0.00244059 0.000814475 0.00025394 0.000433384 6.78098e-05 0.00128319 0.0007984 0.000663954 0.000919029 0.000229612 0.000305659 0.000317644 6.3906e-05
# F6_CAN_LL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0030264 0.00201053 0.00390185 0.00232193 0.00242058 0.0026886 0.00270109 0.00270652 0.00273801 0.00292675 0.00302219 0.00253796 0.00241031 0.00142313 0.0018129 0.00142667 0.00142643 0.00111048 0.0014819 0.00239181 0.00399315 0.000803371
# F7_MOR_LL 0 0 0 0 0 0 0 0 0 0 0 0.00013519 0.000270849 0.000135612 0.000407051 0.000316839 0.000271949 0.000249828 0.00029594 0.0003652 0.000320167 0.000481295 0.000344611 0.000230193 0.000276888 0.000347394 0.000278358 0.000139213 0.000255421 0.00483205 0.00316595 0.00290714 0.00215672 0.00299664 0.00191657 0.00340958 0.00463847 0.00530234 0.00573186 0.00679722 0.000780397 0.0031184 0.00144735 0.00100879 0.000278717 0.00120551 0.0016815 0.0126909 0 0.00203904 0.002118 0.0134963 0.00721324 0.00980401 0.0108149 0.0132829 0.0128723 0.00834556 0.0150379 0.0256426 0.0368642 0.0243818 0.025647 0.0454267 0.0458225 0.0678436 0.0136493
# F8_USA_RR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00703074 0.0113142 0.0068705 0.0175581 0.0571522 0.026188 0.00977093 0.0176366 0.0155462 0.00648531 0.0052918 0.00660729 0.0186935 0.0095213 0.0460647 0.00832228 0.00646027 0.00630869 0.00311239 0.0133092 0.00971502 0.0109817 7.77819e-06 0.0140264 0.0115434 0.0101065 0.00599092 0.00570306 0.0058601 0.00619552 0.00683666 0.00923033 0.00938156 0.00907068 0.00901243 0.00181319
# F9_BEL_LL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000843689 0.00102557 0.00255662 0.00430855 0.00388311 5.09432e-05 2.61223e-05 5.25548e-06
# F10_MOR_PS 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00187552 0.00175479 0.00373886 0.00347619 0.0117696 0.0023679
# F11_CPR_LL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.18955e-05 0 0 0 0 0 0 0.00398827 0.000731364 0.000868079 0.00131432 0.000830383 0.00116865 0.000588059 0.000268967 8.61665e-05 1.73356e-05
# F12_OTH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8.88213e-07 0 0 3.12445e-05 1.91626e-05 1.16709e-05 1.43619e-05 0.00029832 0.000263403 0.000205882 0.000242727 0.000138089 0.00063444 0.00154617 0.000932089 0.00090695 0.0014666 0.00175364 0.00239487 0.00195883 0.00143726 0.00167673 0.00564783 0.00199622 0.00116336 0.00418773 0.00349677 0.00499068 0.00108824 0.00187959 0.0025644 0.0026856 0.00227752 0.00122521 0.000246496
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
        13         1         0         0         0         0  #  S1_USA_LL_Log
        14         1         0         0         0         1  #  S2_USA_LL_Obs
        15         1         0         0         0         0  #  S3_JPN_LL
        16         1         0         0         0         0  #  S4_EU_POR_LL
        17         1         0         0         0         0  #  S5_EU_ESP_LL
        18         1         0         0         0         0  #  S6_CTP_LL
-9999 0 0 0 0 0
#
#_Q_parameters
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -25            25      -6.77573             0             0.05             1         1          0          0          0          0          0          0          0  #  LnQ_base_S1_USA_LL_Log(13)
           -25            25      -6.52246             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_S2_USA_LL_Obs(14)
           -25            25      -8.08215             0             0.05             1         1          0          0          0          0          0          0          0  #  LnQ_base_S3_JPN_LL(15)
           -25            25       -6.2173             0             0.05             1         1          0          0          0          0          0          0          0  #  LnQ_base_S4_EU_POR_LL(16)
           -25            25      -5.90059             0             0.05             1         1          0          0          0          0          0          0          0  #  LnQ_base_S5_EU_ESP_LL(17)
           -25            25      -8.46722             0             0.05             1         1          0          0          0          0          0          0          0  #  LnQ_base_S6_CTP_LL(18)
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
 24 0 0 0 # 1 F1_EU_LL
 24 0 4 0 # 2 F2_JPN_LL
 24 0 3 0 # 3 F3_CTP_LL
 24 0 4 0 # 4 F4_USA_LL
 24 0 4 0 # 5 F5_VEN_LL
 5 0 0 4 # 6 F6_CAN_LL
 5 0 0 1 # 7 F7_MOR_LL
 5 0 0 4 # 8 F8_USA_RR
 5 0 0 5 # 9 F9_BEL_LL
 5 0 0 1 # 10 F10_MOR_PS
 5 0 0 3 # 11 F11_CPR_LL
 5 0 0 3 # 12 F12_OTH
 5 0 0 4 # 13 S1_USA_LL_Log
 5 0 0 4 # 14 S2_USA_LL_Obs
 5 0 0 2 # 15 S3_JPN_LL
 5 0 0 1 # 16 S4_EU_POR_LL
 5 0 0 1 # 17 S5_EU_ESP_LL
 5 0 0 3 # 18 S6_CTP_LL
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
 11 0 0 0 # 1 F1_EU_LL
 11 0 0 0 # 2 F2_JPN_LL
 11 0 0 0 # 3 F3_CTP_LL
 11 0 0 0 # 4 F4_USA_LL
 11 0 0 0 # 5 F5_VEN_LL
 11 0 0 0 # 6 F6_CAN_LL
 11 0 0 0 # 7 F7_MOR_LL
 11 0 0 0 # 8 F8_USA_RR
 11 0 0 0 # 9 F9_BEL_LL
 11 0 0 0 # 10 F10_MOR_PS
 11 0 0 0 # 11 F11_CPR_LL
 11 0 0 0 # 12 F12_OTH
 11 0 0 0 # 13 S1_USA_LL_Log
 11 0 0 0 # 14 S2_USA_LL_Obs
 11 0 0 0 # 15 S3_JPN_LL
 11 0 0 0 # 16 S4_EU_POR_LL
 11 0 0 0 # 17 S5_EU_ESP_LL
 11 0 0 0 # 18 S6_CTP_LL
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   F1_EU_LL LenSelex
          62.5         297.5        138.71        135.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F1_EU_LL(1)
            -6             4      -5.51904            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F1_EU_LL(1)
            -1             9       6.78788           6.7          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F1_EU_LL(1)
            -1             9       7.20099          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F1_EU_LL(1)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F1_EU_LL(1)
            -5             9         -4.96            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F1_EU_LL(1)
# 2   F2_JPN_LL LenSelex
          62.5         297.5       169.941        148.87          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F2_JPN_LL(2)
            -6             4      -4.49702         -4.56          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F2_JPN_LL(2)
            -1             9       7.55227          7.25          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F2_JPN_LL(2)
            -1             9       6.39224          7.61          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F2_JPN_LL(2)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F2_JPN_LL(2)
            -5             9      -2.75508            -5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F2_JPN_LL(2)
           -60           200       -21.842         -2.81            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F2_JPN_LL(2)
           -15            15     -0.783791         -0.93            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F2_JPN_LL(2)
           -15            15       1.26635         -1.15            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F2_JPN_LL(2)
           -15            15      -2.15177             0            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F2_JPN_LL(2)
           -15            15      0.780966          0.62            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F2_JPN_LL(2)
# 3   F3_CTP_LL LenSelex
          62.5         297.5       166.352        159.98          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F3_CTP_LL(3)
            -6             4       -2.1352            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F3_CTP_LL(3)
            -1             9       6.73051          6.81          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F3_CTP_LL(3)
            -1             9       6.80631          7.08          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F3_CTP_LL(3)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F3_CTP_LL(3)
            -5             9      -4.03119            -5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F3_CTP_LL(3)
          -200           200      -25.6791         -2.81            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Peak_F3_CTP_LL(3)
           -15            15     -0.876797         -0.93            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Ascend_F3_CTP_LL(3)
           -15            15      0.186427         -1.15            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Descend_F3_CTP_LL(3)
           -15            15       1.12049             0            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Final_F3_CTP_LL(3)
           -15            15      0.552691          0.62            50             0          5          0          0          0          0          0          0          0  #  SzSel_Male_Scale_F3_CTP_LL(3)
# 4   F4_USA_LL LenSelex
          62.5         297.5       165.505        127.99          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F4_USA_LL(4)
            -6             4      -2.88348         -5.84          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F4_USA_LL(4)
            -1             9       7.41005          7.33          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F4_USA_LL(4)
            -1             9       6.51128          8.08          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F4_USA_LL(4)
            -5             9      -1.36928          -2.5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F4_USA_LL(4)
            -5             9       -2.4305            -5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F4_USA_LL(4)
           -40           200       -25.243         -2.81            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F4_USA_LL(4)
           -15            15      -1.53889         -0.93            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F4_USA_LL(4)
           -15            15       1.38621         -1.15            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F4_USA_LL(4)
           -15            15      -2.44013             0            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F4_USA_LL(4)
           -15            15      0.824405          0.62            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F4_USA_LL(4)
# 5   F5_VEN_LL LenSelex
          62.5         297.5       177.287        167.54          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F5_VEN_LL(5)
            -6             4       -5.1778            -6          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F5_VEN_LL(5)
            -1             9       7.24856          6.81          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F5_VEN_LL(5)
            -1             9       6.78851          7.08          0.05             1          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F5_VEN_LL(5)
            -5             9            -5            -5          0.05             1         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F5_VEN_LL(5)
            -5             9      -2.19589            -5          0.05             1          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F5_VEN_LL(5)
          -200           200        20.196         -2.81            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F5_VEN_LL(5)
           -15            15       0.94367         -0.93            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F5_VEN_LL(5)
           -15            15      0.201728         -1.15            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F5_VEN_LL(5)
           -15            15      -1.83318             0            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F5_VEN_LL(5)
           -15            15      0.724023          0.62            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F5_VEN_LL(5)
# 6   F6_CAN_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F6_CAN_LL(6)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F6_CAN_LL(6)
# 7   F7_MOR_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F7_MOR_LL(7)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F7_MOR_LL(7)
# 8   F8_USA_RR LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F8_USA_RR(8)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F8_USA_RR(8)
# 9   F9_BEL_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F9_BEL_LL(9)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F9_BEL_LL(9)
# 10   F10_MOR_PS LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F10_MOR_PS(10)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F10_MOR_PS(10)
# 11   F11_CPR_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F11_CPR_LL(11)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F11_CPR_LL(11)
# 12   F12_OTH LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F12_OTH(12)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F12_OTH(12)
# 13   S1_USA_LL_Log LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S1_USA_LL_Log(13)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S1_USA_LL_Log(13)
# 14   S2_USA_LL_Obs LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S2_USA_LL_Obs(14)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S2_USA_LL_Obs(14)
# 15   S3_JPN_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S3_JPN_LL(15)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S3_JPN_LL(15)
# 16   S4_EU_POR_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S4_EU_POR_LL(16)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S4_EU_POR_LL(16)
# 17   S5_EU_ESP_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S5_EU_ESP_LL(17)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S5_EU_ESP_LL(17)
# 18   S6_CTP_LL LenSelex
             0            10             1             1            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S6_CTP_LL(18)
            10           100            66            66            25             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S6_CTP_LL(18)
# 1   F1_EU_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F1_EU_LL(1)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F1_EU_LL(1)
# 2   F2_JPN_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F2_JPN_LL(2)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F2_JPN_LL(2)
# 3   F3_CTP_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F3_CTP_LL(3)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F3_CTP_LL(3)
# 4   F4_USA_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F4_USA_LL(4)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F4_USA_LL(4)
# 5   F5_VEN_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F5_VEN_LL(5)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F5_VEN_LL(5)
# 6   F6_CAN_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F6_CAN_LL(6)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F6_CAN_LL(6)
# 7   F7_MOR_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F7_MOR_LL(7)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F7_MOR_LL(7)
# 8   F8_USA_RR AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F8_USA_RR(8)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F8_USA_RR(8)
# 9   F9_BEL_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F9_BEL_LL(9)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F9_BEL_LL(9)
# 10   F10_MOR_PS AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F10_MOR_PS(10)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F10_MOR_PS(10)
# 11   F11_CPR_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F11_CPR_LL(11)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F11_CPR_LL(11)
# 12   F12_OTH AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F12_OTH(12)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F12_OTH(12)
# 13   S1_USA_LL_Log AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S1_USA_LL_Log(13)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S1_USA_LL_Log(13)
# 14   S2_USA_LL_Obs AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S2_USA_LL_Obs(14)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S2_USA_LL_Obs(14)
# 15   S3_JPN_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S3_JPN_LL(15)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S3_JPN_LL(15)
# 16   S4_EU_POR_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S4_EU_POR_LL(16)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S4_EU_POR_LL(16)
# 17   S5_EU_ESP_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S5_EU_ESP_LL(17)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S5_EU_ESP_LL(17)
# 18   S6_CTP_LL AgeSelex
             0            10             0             0            25             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S6_CTP_LL(18)
            10           100            30             0            25             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S6_CTP_LL(18)
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
      4      1    0.0483
      5      1         0
      6      1         0
      4      2    0.0569
      5      2         0
      6      2         0
      4      3    0.0396
      5      3         0
      6      3         0
      4      4    0.0997
      5      4         0
      6      4         0
      4      5    0.2537
      5      5         0
      6      5         0
      5      6         0
      6      6         0
      5      7         0
      6      7         0
      5      8         0
      6      8         0
      5      9         0
      6      9         0
      5     10         0
      6     10         0
      5     11         0
      6     11         0
      5     12         0
      6     12         0
      5     13         0
      6     13         0
      5     14         0
      6     14         0
      1     15    0.1459
      5     15         0
      6     15         0
      1     16    0.0578
      5     16         0
      6     16         0
      1     17    0.0886
      5     17         0
      6     17         0
      1     18     0.251
      5     18         0
      6     18         0
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 19 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 1 13 1 1 1
 1 14 1 0 1
 1 15 1 1 1
 1 16 1 1 1
 1 17 1 1 1
 1 18 1 1 1
 4 1 1 1 1
 4 2 1 1 1
 4 3 1 1 1
 4 4 1 1 1
 4 5 1 1 1
 9 1 1 1 1
 9 2 1 1 1
 9 7 1 1 1
 10 1 1 1 1
 11 1 1 1 1
 12 1 1 1 1
 13 1 1 1 1
 17 1 1 0 1
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
#  0 #_CPUE/survey:_12
#  1 #_CPUE/survey:_13
#  0 #_CPUE/survey:_14
#  1 #_CPUE/survey:_15
#  1 #_CPUE/survey:_16
#  1 #_CPUE/survey:_17
#  1 #_CPUE/survey:_18
#  1 #_lencomp:_1
#  1 #_lencomp:_2
#  1 #_lencomp:_3
#  1 #_lencomp:_4
#  1 #_lencomp:_5
#  0 #_lencomp:_6
#  0 #_lencomp:_7
#  0 #_lencomp:_8
#  0 #_lencomp:_9
#  0 #_lencomp:_10
#  0 #_lencomp:_11
#  0 #_lencomp:_12
#  0 #_lencomp:_13
#  0 #_lencomp:_14
#  0 #_lencomp:_15
#  0 #_lencomp:_16
#  0 #_lencomp:_17
#  0 #_lencomp:_18
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

