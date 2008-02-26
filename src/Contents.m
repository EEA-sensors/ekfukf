% EKF/UKF toolbox for Matlab 7.x
% Version 1.2, February 8. 2008
%
% Copyright (C) 2005-2008 Simo Särkkä, <simo.sarkka@hut.fi>
%               2007-2008 Jouni Hartikainen <jmjharti@cc.hut.fi>
% History:      
%   04.09.2007 JH & SS Updated for version 1.1
%   06.08.2007 JH Updated for version 1.0
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.
% 
%
% Kalman filtering
%   KF_PREDICT    Perform Kalman Filter prediction step
%   KF_UPDATE     Kalman Filter update step
%   KF_LHOOD      Kalman Filter measurement likelihood
%   RTS_SMOOTH    Rauch-Tung-Striebel Smoother
%   TF_SMOOTH     Smoother based on combination of two Kalman filters
%
% Extended Kalman filtering
%   EKF_PREDICT1  1st order Extended Kalman Filter prediction step
%   EKF_UPDATE1   1st order Extended Kalman Filter update step
%   EKF_PREDICT2  2nd order Extended Kalman Filter prediction step
%   EKF_UPDATE2   2nd order Extended Kalman Filter update step
%   ERTS_SMOOTH1  1st order Extended RTS Smoother
%   ETF_SMOOTH1   Smoother based on two 1. order extended Kalman filters           
%
% Unscented transform / Unscented Kalman filtering
%   UT_WEIGHTS    Generate weights for sigma points using the summation form
%   UT_MWEIGTS    Generate weights for sigma points using the matrix form
%   UT_SIGMAS     Generate Sigma Points for Unscented Transformation
%   UT_TRANSFORM  Makes the Unscented Transformation of x and y
%   UKF_PREDICT1  Nonaugmented UKF prediction step
%   UKF_UPDATE1   Nonaugmented UKF update step
%   UKF_PREDICT2  Augmented (state and process noise) UKF prediction step 
%   UKF_UPDATE2   Augmented (state and measurement noise) UKF update step 
%   UKF_PREDICT3  Augmented (state, process and measurement noise) UKF prediction step
%   UKF_UPDATE3   Augmented (state, process and measurement noise) UKF update step
%   URTS_SMOOTH1  Nonaugmented unscented RTS-smoother
%   URTS_SMOOTH2  Augmented unscented RTS-smoother
%   UTF_SMOOTH    Smoother based on combination of two unscented Kalman filters
%
% Multiple Model Filtering
%   IMM_PREDICT   IMM filter prediction step
%   IMM_UPDATE    IMM filter update step
%   IMM_SMOOTH    IMM smoothing
%   EIMM_PREDICT  IMM-EKF filter prediction step
%   EIMM_UPDATE   IMM-EKF filter update step
%   EIMM_SMOOTH   IMM-EKF smoothing
%   UIMM_PREDICT  IMM-UKF filter prediction step
%   UIMM_UPDATE   IMM-UKF filter update step
%   UIMM_SMOOTH   IMM-UKF smoothing
%
%
% Misc.
%   GAUSS_PDF     Multivariate Gaussian PDF
%   GAUSS_RND     Multivariate Gaussian random variables
%   LTI_INT       Integrate LTI ODE with Gaussian Noise
%   LTI_DISC      Discretize LTI ODE with Gaussian Noise
%   RK4           Runge-Kutta integration
%   DER_CHECK     Check derivatives using finite differences
%   SCHOL         Positive semidefinite matrix Cholesky factorization
%   RESAMPSTR     Stratified resampling
%
% /DEMOS/ 
%
%   /KF_CWPA_DEMO/             
%      KF_CWPA_DEMO       CWPA model demonstration with Kalman filter
%
%   /EKF_SINE_DEMO/          
%      EKF_SINE_F         Dynamic model function (needed by the augmented UKF)
%      EKF_SINE_H         Measurement model function
%      EKF_SINE_DH_DX     1st order derivative of the measurement model
%      EKF_SINE_D2H_DX2   2nd order derivative of the measurement model
%      EKF_SINE_DEMO      Random Sine Signal demonstration
%
%   /UNGM_DEMO/           
%      UNGM_F             Dynamic model function
%      UNGM_DF_DX         1st order derivative of the dynamic model
%      UNGM_D2F_DX2       2nd order derivative of the dynamic model (not used)
%      UNGM_H             Measurement model function
%      UNGM_DH_DX         1st order derivative of the measurement model
%      UNGM_D2H_DX2       2nd order derivative of the measurement model (not used)
%      UNGM_DEMO          UNGM model demonstration
%
%   /BOT_DEMO/            
%      BOT_H              Measurement model function
%      BOT_DH_DX          1st order derivative of the measurement model 
%      BOT_D2H_DX2        2nd order derivative of the measurement model  
%      BOT_DEMO_ALL       BOT demo with EKF and UKF
%      EKFS_BOT_DEMO      BOT demo with EKF
%      UKFS_BOT_DEMO      BOT demo with UKF
%
%   /REENTRY_DEMO/        
%      REENTRY_F          Dynamic model function
%      REENTRY_DF         Derivative of the dynamic model
%      REENTRY_H          Measurement model function
%      REENTRY_DH         Derivative of the measurement model
%      REENTRY_IF         Inverse prediction of the dynamic model
%      REENTRY_COND       Generates condition numbers for simulation data
%      MAKE_REENTRY_DATA  Generates the simulation data for reentry dynamics 
%      REENTRY_DEMO       Reentry Vehicle Tracking demonstration
%  
%   /IMM_DEMO/
%      IMM_DEMO           Tracking a Target with Simple Manouvers demonstration
%
%   /EIMM_DEMO/
%      F_TURN             Dynamic model function for the coordinated turn model
%      F_TURN_DX          Jacobian of the coordinated turn model's dynamic model
%      F_TURN_INV         Inverse dynamics of the coordinated turn model
%      CT_DEMO            Coordinated Turn Model demonstration
%      BOT_H              Measurement model function
%      BOT_DH_DX          1st order derivative of the measurement model 
%      BOT_D2H_DX2        2nd order derivative of the measurement model  
%      BOTM_DEMO          Bearings Only Tracking of a Manouvering Target Demonstration
%
% Demos currently included in the toolbox, but not documented:
%
% /KF_SINE_DEMO/           
%      KF_SINE_DEMO       Sine signal demonstration with Kalman filter
