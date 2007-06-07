% UKF/EKF toolbox for Matlab 7.x
% Version 0.1, February 19 2007
%
% Copyright (C) 2005 Simo Särkkä, <simo.sarkka@hut.fi>
%               2007 Jouni Hartikainen <jmjharti@cc.hut.fi>
% History:      
%   07.06.2007 JH Initial version for the toolbox. Modified from the SS's
%                 previous work.
%         
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
%   FBF_SMOOTH    Forward-Backward Smoother
%
% Extended Kalman filtering
%   EKF_PREDICT1  1st order Extended Kalman Filter prediction step
%   EKF_UPDATE1   1st order Extended Kalman Filter update step
%   EKF_PREDICT2  2nd order Extended Kalman Filter prediction step
%   EKF_UPDATE2   2nd order Extended Kalman Filter update step
%   ERTS_SMOOTH1  1st order Extended RTS Smoother
%   EFBF_SMOOTH1  1st order Forward-Backward Smoother            
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
%   UKF_PREDICT3  Augmented (state,process and measurement noise) UKF prediction step
%   UKF_UPDATE3   Augmented (state,process and measurement noise) UKF update step
%
% Misc.
%   GAUSS_PDF  Multivariate Gaussian PDF
%   GAUSS_RND  Multivariate Gaussian random variables
%   LTI_INT    Integrate LTI ODE with Gaussian Noise
%   LTI_DISC   Discretize LTI ODE with Gaussian Noise
%   RK4        Runge-Kutta integration
%   DER_CHECK  Check derivatives using finite differences
%   SCHOL      Positive semidefinite matrix Cholesky factorization
%
% Demos 
%   BOT_DEMO - Bearings Only Tracking demonstration
%      AZ_H                Measurement model function
%      AZ_DH_DX            1st order derivative of the measurement model 
%      AZ_D2H_DX2          2nd order derivative of the measurement model  
%      BOT_DEMO_ALL        BOT demo with EKF and UKF
%      EKFS_BOT_DEMO       BOT demo with EKF
%      UKFS_BOT_DEMO       BOT demo with UKF
%
%   DREENTRY - Reentry Vehicle Tracking demonstration
%      DREENTRY_A          Dynamic model function
%      DREENTRY_DA         Derivative of the dynamic model
%      DREENTRY_H          Measurement model function
%      DREENTRY_DH         Derivative of the measurement model
%      DREENTRY_IA         Inverse prediction of the dynamic model
%      DREENTRY_COND       Generates condition numbers for simulation data
%      MAKE_DREENTRY_DATA  Generates the simulation data for reentry dynamics 
%      DREENTRY_DEMO       Reentry demonstration
%
%   EKF_SINE_DEMO - Random Sine Signal Tracking demonstration   
%      EKF_SINE_F          Dynamic model function (needed by the augmented UKF)
%      EKF_SINE_H          Measurement model function
%      EKF_SINE_DH_DX      1st order derivative of the measurement model
%      EKF_SINE_D2H_DX2    2nd order derivative of the measurement model
%      EKF_SINE_DEMO       Random Sine Signal demonstration
%
%   KF_CWPA_DEMO  - CWPA model demonstration with Kalman filter
%      KF_CWPA_DEMO        CWPA model demonstration
%
%   KF_SINE_DEMO  - Sine Signal demonstration with Kalman filter   
%      KF_SINE_DEMO        Sine signal demonstration
%
%   UNGM_DEMO     - UNGM model demonstration
%      UNGM_F              Dynamic model function
%      UNGM_DF_DX          1st order derivative of the dynamic model
%      UNGM_D2F_DX2        2nd order derivative of the dynamic model (not used)
%      UNGM_H              Measurement model function
%      UNGM_DH_DX          1st order derivative of the measurement model
%      UNGM_D2H_DX2        2nd order derivative of the measurement model (not used)
%      UNGM_DEMO           UNGM model demonstration