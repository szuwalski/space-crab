//
//  spatial_size_1.cpp
//  Created by Jie Cao on 10/27/17 (modified from MIST developed by Jim Thorson).
//

#include <TMB.hpp>

// Posfun
template<class Type>
Type posfun(Type x, Type lowerlimit, Type &pen){
    pen += CppAD::CondExpLt(x,lowerlimit,Type(0.01)*pow(x-lowerlimit,2),Type(0));
    return CppAD::CondExpGe(x,lowerlimit,x,lowerlimit/(Type(2)-x/lowerlimit));
}

// 2nd power of a number
template<class Type>
Type square(Type x){ return x*x; }

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
    return R_IsNA(asDouble(x));
}

// function for logistic transform
template<class Type>
Type plogis(Type x){
    return 1.0 / (1.0 + exp(-x));
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
    Type Return;
    if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
    if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
    return Return;
}

// dzinflognorm
template<class Type>
Type dzinflognorm(Type x, Type meanlog, Type encounter_prob, Type log_notencounter_prob, Type sdlog, int give_log=false){
    Type Return;
    if(x==0){
        if(give_log==false) Return = 1.0 - encounter_prob;
        if(give_log==true){
            if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
            if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
        }
    }else{
        if(give_log==false) Return = encounter_prob * dlognorm( x, meanlog, sdlog, false );
        if(give_log==true) Return = log(encounter_prob) + dlognorm( x, meanlog, sdlog, true );
    }
    return Return;
}

// dzinfgamma, shape = 1/CV^2, scale = mean*CV^2
template<class Type>
Type dzinfgamma(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
    Type Return;
    if(x==0){
        if(give_log==false) Return = 1.0 - encounter_prob;
        if(give_log==true){
            if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
            if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
        }
    }else{
        if(give_log==false) Return = encounter_prob * dgamma( x, pow(cv,-2), posmean*pow(cv,2), false );
        if(give_log==true) Return = log(encounter_prob) + dgamma( x, pow(cv,-2), posmean*pow(cv,2), true );
    }
    return Return;
}

// dzinfnorm
template<class Type>
Type dzinfnorm(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
    Type Return;
    if(x==0){
        if(give_log==false) Return = 1.0 - encounter_prob;
        if(give_log==true){
            if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
            if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
        }
    }else{
        if(give_log==false) Return = encounter_prob * dnorm( x, posmean, posmean*cv, false );
        if(give_log==true) Return = log(encounter_prob) + dnorm( x, posmean, posmean*cv, true );
    }
    return Return;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    
    using namespace R_inla;
    using namespace Eigen;
    using namespace density;
    
    DATA_IVECTOR( Options_vec );
    // Slot 0 -- option for diagonal covariance (i.e., independence among species)
    // Slot 1 -- method for spatial variation; 0=SPDE; 1=AR1 precision matrix
    DATA_IVECTOR( ObsModel_p );      // Slot 1-n_p -- distribution of data: 0=Poisson; 1=Lognormal; 2=Zero-inflated lognormal; 3=lognormal-Poisson; 4=Normal
    
    // Indices
    DATA_INTEGER(n_i);       // Total number of observations (i)
    //DATA_INTEGER(n_s);       // Number of stations (s)
    DATA_INTEGER(n_t);       // Number of years (t)
    DATA_INTEGER(n_k);       // Number of knots (k)
    DATA_INTEGER(n_p);       // Number of size classes (p)
    DATA_INTEGER(n_j);   // Number of dynamic factors in process error (j)
    
    // Data
    DATA_VECTOR( b_i );      // Count for observation
    DATA_FACTOR( p_i );      // Size class for observation * 2
    DATA_FACTOR( s_i );      // Site for observation
    DATA_FACTOR( t_i );      // Year for observation
    DATA_VECTOR( R_sex );    // sex ratio for recruitment first slot -male; second slot -female
    DATA_VECTOR( R_size );   // proportion of recruitment for each size bin
    DATA_SCALAR( M_frac);    // fraction of year for pulse fishery
    DATA_VECTOR( immature);  // 1 - maturity at size
    DATA_VECTOR( mature);    // maturity at size
    
    DATA_ARRAY( M_pkt );     // natural mortality
    DATA_MATRIX( growthmale_tran );  //growth transition matrix for males
    DATA_MATRIX( growthfemale_tran ); //growth transition matrix for females
    
    // Harvest
    DATA_ARRAY( c_pkt );           // Harvest for each location, year, and species
    DATA_SCALAR( catchlogsd );     // Logsd for harvest
    
    DATA_VECTOR( a_k );        // Area for each "real" stratum(km^2) in each stratum (zero for all knots k not associated with stations s)
    //DATA_MATRIX( Z_kl );        // Derived quantity matrix
    
    // Aniso objects
    DATA_STRUCT(spde, spde_aniso_t);
    
    // Sparse matrices for precision matrix of 2D AR1 process
    // Q = M0*(1+rho^2)^2 + M1*(1+rho^2)*(-rho) + M2*rho^2
    DATA_SPARSE_MATRIX(M0);
    DATA_SPARSE_MATRIX(M1);
    DATA_SPARSE_MATRIX(M2);
    
    // Fixed effects
    PARAMETER_VECTOR(Hinput_z); // Anisotropy parameters
    PARAMETER_VECTOR(logkappa_z);         // Controls range of spatial variation.  First n_p slots are independent for each spatial component (but can be fixed to be equal).  slot n_p is for spatio-temporal (SDFA) component
    PARAMETER_VECTOR(R_mu);   // Mean R
    PARAMETER_VECTOR(phi_p);              // equilibrium N for each size class, same for both sex
    PARAMETER_VECTOR(logMargSigmaR);        // log-inverse SD of Alpha  // logtauR
    PARAMETER_VECTOR(L_val);    // Values in loadings matrix
    PARAMETER_MATRIX(logsigma_pz);
    
    // Harvest rates
    PARAMETER_ARRAY( logF_male_pkt );  // log-instantaneous mortality rate
    PARAMETER_ARRAY( logF_female_pkt); // log-instantaneous mortality rate
    
    // Random effects
    PARAMETER_ARRAY(d_pkt);  // Spatial process variation - density
    PARAMETER_MATRIX(Rinput_kt);  // Spatial process variation - recruitment
    PARAMETER_VECTOR(delta_i);
    
    // global stuff
    //int n_l = Z_kl.row(0).size();
    Type jnll = 0;
    vector<Type> jnll_comp(6);
    jnll_comp.setZero();
    matrix<Type> Identity_pp(n_p*2, n_p*2);
    Identity_pp.setIdentity();
    Type pos_penalty =0;
    
    // Anisotropy elements
    matrix<Type> H(2,2);
    H(0,0) = exp(Hinput_z(0));
    H(1,0) = Hinput_z(1);
    H(0,1) = Hinput_z(1);
    H(1,1) = (1+Hinput_z(1)*Hinput_z(1)) / exp(Hinput_z(0));
    
    // Covariance for process error
    matrix<Type> L_pj(n_p*2, n_j);
    matrix<Type> L_jp(n_j, n_p*2);
    L_pj.setZero();
    matrix<Type> Cov_pp(n_p*2, n_p*2);
    Cov_pp.setZero();
    int Count = 0;
    // Trimmed-Cholesky "factor model"
    if( Options_vec(0)==0 ){
        // Assemble the loadings matrix (lower-triangular, loop through rows then columns)
        for(int j=0; j<n_j; j++){
            for(int p=0; p<n_p*2; p++){
                if(j<=p){
                    L_pj(p,j) = L_val(Count);
                    Count++;
                }
            }}
        // Calculate the covariance
        L_jp = L_pj.transpose();
        Cov_pp = L_pj*L_jp + Type(0.000001)*Identity_pp;  // additive constant to make Cov_pp invertible
    }
    // Independent among size classes
    if( Options_vec(0)==1 ){
        for(int p=0; p<n_p*2; p++){
            Cov_pp(p,p) = exp( L_val(p) );
        }
    }
    
    // Derived quantities related to spatial variation (logtauR=Spatiotemporal1;  logtauE=Spatiotemporal2)
    // calculation of logtauE_p and Range_pz depends upon whether we're treating size classes as independent or not
    vector<Type> Range_pz(n_p+1); //fist slot for recruitment, the rest slots for spatiotemporal2
    vector<Type> logtauE_p(n_p);
    vector<Type> logtauR(1);
    vector<Type> MargSigmaR(1);
    
        logtauR(0) = 0.5*log(4*M_PI) - logMargSigmaR(0) - logkappa_z(n_p);
        Range_pz(0) = sqrt(8) / exp( logkappa_z(n_p) );
        MargSigmaR(0) = exp( logMargSigmaR(0) );
    
    for(int p=0; p<n_p; p++){
        // Not independent
        if( Options_vec(0)==0 ){
            logtauE_p(p) = log(1/( exp(logkappa_z(n_p)) * sqrt(4*M_PI)) );
            Range_pz(p+1) = sqrt(8) / exp( logkappa_z(n_p) );
        }
        // Independent among species
        if( Options_vec(0)==1 ){
            logtauE_p(p) = log(1/( exp(logkappa_z(p)) * sqrt(4*M_PI)) );
            Range_pz(p+1) = sqrt(8) / exp( logkappa_z(p) );
        }
    }
    
    // Derived quantities related to GMRF
    Eigen::SparseMatrix<Type> Q_spatiotemporal1;
    Eigen::SparseMatrix<Type> Q_spatiotemporal2;
    MVNORM_t<Type> nll_mvnorm(Cov_pp);
    REPORT( Q_spatiotemporal1 );
    REPORT( Q_spatiotemporal2 );
    
    // Transform random fields
    matrix<Type> R_kt(n_k, n_t);
    for(int k=0; k<n_k; k++){
        for(int t=0; t<n_t; t++){
            R_kt(k,t) = Rinput_kt(k,t)/exp(logtauR(0)) + R_mu(t);
        }}
    
    //Calculate predicted density in year t+1 given model and density in year t
    vector<Type> tmp1_p(n_p);
    vector<Type> tmp1a_p(n_p);
    vector<Type> tmp2_p(n_p);
    vector<Type> tmp3_p(n_p);
    vector<Type> tmp3a_p(n_p);
    vector<Type> tmp4_p(n_p);
    array<Type> d1_pkt(n_p, n_k, n_t); //
    array<Type> d2_pkt(n_p, n_k, n_t); //
    array<Type> dhat1_pkt(n_p, n_k, n_t);       // male density
    array<Type> dhat2_pkt(n_p, n_k, n_t);       // female density
    array<Type> dpred1_pkt(n_p, n_k, n_t);      // Predicted male density for each year
    array<Type> dpred2_pkt(n_p, n_k, n_t);      // Predicted female density for each year
    array<Type> dpred_pkt(2*n_p, n_k, n_t);      // rbind(dpred1,dpred2) for each year
    array<Type> F_male_pkt(n_p, n_k, n_t);        // F
    array<Type> F_female_pkt(n_p, n_k, n_t);        // F
    array<Type> chat_pkt( n_p, n_k, n_t );      // Predicted catch-at-size for each year
    chat_pkt.setZero();
    
    // First year
    for(int k=0; k<n_k; k++){
        for(int p=0; p<n_p; p++){
            dpred1_pkt(p,k,0) = phi_p(p) + R_kt(k,0)*R_sex(0)*R_size(p); /////// male - 1; female -2
            dpred2_pkt(p,k,0) = phi_p(p) + R_kt(k,0)*R_sex(1)*R_size(p);
        }}
    
    for(int t=0; t<n_t; t++){
    for(int k=0; k<n_k; k++){
        
        for (int p=0; p<n_p; p++){
            d1_pkt(p,k,t) = d_pkt(p,k,t);
            d2_pkt(p,k,t) = d_pkt(n_p+p,k,t);
        }
    
        //F_male_pkt.col(t).col(k) = exp( logF_male_pkt.col(t).col(k) );
        //F_female_pkt.col(t).col(k) = exp( logF_female_pkt.col(t).col(k) );
        
        //tmp1_p = growthmale_tran * (d1_pkt.col(t).col(k) * immature * exp( -M_pkt.col(t).col(k)-F_male_pkt.col(t).col(k) )); // immature N
        //tmp2_p = d1_pkt.col(t).col(k) * mature * exp( -M_pkt.col(t).col(k)-F_male_pkt.col(t).col(k) ); // mature N
        //dhat1_pkt.col(t).col(k) = tmp1_p + tmp2_p;
    
        //tmp3_p = growthfemale_tran * (d2_pkt.col(t).col(k) * immature * exp( -M_pkt.col(t).col(k)-F_female_pkt.col(t).col(k) )); // immature N
        //tmp4_p = d2_pkt.col(t).col(k) * mature * exp( -M_pkt.col(t).col(k)-F_female_pkt.col(t).col(k) ); // mature N
        //dhat2_pkt.col(t).col(k) = tmp3_p + tmp4_p;
        
        //chat_pkt.col(t).col(k) = a_k(k) * (1 - exp( -F_male_pkt.col(t).col(k) )) * d1_pkt.col(t).col(k) * exp( -M_frac * M_pkt.col(t).col(k) ); // catch
        
        for (int p=0; p<n_p; p++){
            F_male_pkt(p,k,t) = exp( logF_male_pkt(p,k,t) );
            F_female_pkt(p,k,t) = exp( logF_female_pkt(p,k,t) );
            tmp1a_p(p) = d1_pkt(p,k,t) * exp(-M_pkt(p,k,t)-F_male_pkt(p,k,t)) * immature(p);
            tmp2_p(p) = d1_pkt(p,k,t) * exp(-M_pkt(p,k,t)-F_male_pkt(p,k,t)) * mature(p);
            tmp3a_p(p) = d2_pkt(p,k,t) * exp(-M_pkt(p,k,t)-F_female_pkt(p,k,t)) * immature(p);
            tmp4_p(p) = d2_pkt(p,k,t) * exp(-M_pkt(p,k,t)-F_female_pkt(p,k,t)) * mature(p);
        }
        
        tmp1_p = growthmale_tran * tmp1a_p;
        tmp3_p = growthfemale_tran * tmp3a_p;
        
        for (int p=0; p<n_p; p++){
            dhat1_pkt(p,k,t) = posfun((tmp1_p(p) + tmp2_p(p)),Type(1e-10),pos_penalty);
            dhat2_pkt(p,k,t) = posfun((tmp3_p(p) + tmp4_p(p)),Type(1e-10),pos_penalty);
            chat_pkt(p,k,t)  = a_k(k) * (1 - exp( -F_male_pkt(p,k,t) )) * d1_pkt(p,k,t) * exp( -M_frac * M_pkt(p,k,t) );
        }
        
        if( (t+1)<n_t ){
            dpred1_pkt.col(t+1).col(k) = dhat1_pkt.col(t).col(k) + R_kt(k,t+1)*R_sex(0)*R_size;
            dpred2_pkt.col(t+1).col(k) = dhat2_pkt.col(t).col(k) + R_kt(k,t+1)*R_sex(1)*R_size;
        }
        
        for (int p=0; p<n_p; p++){
            dpred_pkt(p,k,t) = dpred1_pkt(p,k,t);
            dpred_pkt(n_p+p,k,t) = dpred2_pkt(p,k,t);
        }
    }}
    REPORT( pos_penalty )
    
    jnll_comp(4) += 1000 * pos_penalty;
    
    // Probability of observed harvest, or of logF if subtracting harvest using posfun to ensure positive biomass
    for(int k=0; k<n_k; k++){
        for(int p=0; p<n_p; p++){
            for(int t=0; t<n_t; t++){
                if( c_pkt(p,k,t)>0 ){
                    jnll_comp(0) -= dnorm( log(c_pkt(p,k,t)), log(chat_pkt(p,k,t)), catchlogsd, true );
                }
            }}}

    // Probability of random fields
    array<Type> Epsilon_kp(n_k, 2*n_p);
    
    // R (spatio-temporal variation in recruitment)
    for(int t=0; t<n_t; t++){
        if( Options_vec(1)==0 ) Q_spatiotemporal1 = Q_spde(spde, exp(logkappa_z(n_p)), H);
        // FOR EXPLORATION SEE: C:\Users\James.Thorson\Desktop\Project_git\2016_Spatio-temporal_models\Week 6 -- 2D spatial models\Lecture\autoregressive_grid_V1.cpp
        if( Options_vec(1)==1 ) Q_spatiotemporal1 = exp(2*logtauR(0)) * (M0*pow(1+exp(logkappa_z(n_p)*2),2) + M1*(1+exp(logkappa_z(n_p)*2))*(-exp(logkappa_z(n_p))) + M2*exp(logkappa_z(n_p)*2));
        jnll_comp(1) += GMRF(Q_spatiotemporal1)(Rinput_kt.col(t));
    }
        
    // Epsilon (spatio-temporal process errors)
    for(int t=0; t<n_t; t++){
        for(int k=0; k<n_k; k++){
            for(int p=0; p<2*n_p; p++){
                Epsilon_kp(k,p) = d_pkt(p,k,t) - dpred_pkt(p,k,t);
            }}
        
        // Depends upon whether using independence or not
        if( Options_vec(0)==0 ){
            if( Options_vec(1)==0 ) Q_spatiotemporal2 = Q_spde(spde, exp(logkappa_z(n_p)), H);
            if( Options_vec(1)==1 ) Q_spatiotemporal2 = exp(2*logtauE_p(0)) * (M0*pow(1+exp(logkappa_z(n_p)*2),2) + M1*(1+exp(logkappa_z(n_p)*2))*(-exp(logkappa_z(n_p))) + M2*exp(logkappa_z(n_p)*2));
            jnll_comp(2) += SCALE(SEPARABLE(nll_mvnorm, GMRF(Q_spatiotemporal2)), exp(-logtauE_p(0)))(Epsilon_kp);
        }
        if( Options_vec(0)==1 ){
            for(int p=0; p<n_p; p++){
                if( Options_vec(1)==0 ) Q_spatiotemporal2 = Q_spde(spde, exp(logkappa_z(p)), H);
                if( Options_vec(1)==1 ) Q_spatiotemporal2 = exp(2*logtauE_p(p)) * (M0*pow(1+exp(logkappa_z(p)*2),2) + M1*(1+exp(logkappa_z(p)*2))*(-exp(logkappa_z(p))) + M2*exp(logkappa_z(p)*2));
                jnll_comp(2) += SCALE( GMRF(Q_spatiotemporal2), exp(L_val(p) - logtauE_p(p)))(Epsilon_kp.col(p));
            }
            
            for(int p=n_p; p<2*n_p; p++){
                if( Options_vec(1)==0 ) Q_spatiotemporal2 = Q_spde(spde, exp(logkappa_z(p-n_p)), H);
                if( Options_vec(1)==1 ) Q_spatiotemporal2 = exp(2*logtauE_p(p-n_p)) * (M0*pow(1+exp(logkappa_z(p-n_p)*2),2) + M1*(1+exp(logkappa_z(p-n_p)*2))*(-exp(logkappa_z(p-n_p))) + M2*exp(logkappa_z(p-n_p)*2));
                jnll_comp(2) += SCALE( GMRF(Q_spatiotemporal2), exp(L_val(p) - logtauE_p(p-n_p)))(Epsilon_kp.col(p));
            }
        }
    }
   
    // Probability of observations
    vector<Type> logchat_i(n_i);
    vector<Type> jnll_i(n_i);
    jnll_i.setZero();
    Type encounterprob;
    Type log_notencounterprob;
        
    for(int i=0; i<n_i; i++){
        logchat_i(i) = log(d_pkt(p_i(i),s_i(i),t_i(i)));
        if( !isNA(b_i(i)) ){
            if( ObsModel_p(p_i(i))==0 ) jnll_i(i) = -1 * dpois( b_i(i), exp(logchat_i(i)), true );
            if( ObsModel_p(p_i(i))==1 ) jnll_i(i) = -1 * dlognorm( b_i(i), logchat_i(i)-pow(exp(logsigma_pz(p_i(i),0)),2)/2, exp(logsigma_pz(p_i(i),0)), true );
            if( ObsModel_p(p_i(i))==2 ){
                encounterprob = ( 1.0 - exp(-1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1))) );
                log_notencounterprob = -1 * exp(logchat_i(i)) * exp(logsigma_pz(p_i(i),1));
                jnll_i(i) = -1 * dzinflognorm( b_i(i), logchat_i(i)-log(encounterprob), encounterprob, log_notencounterprob, exp(logsigma_pz(p_i(i),0)), true);
                }
            if( ObsModel_p(p_i(i))==3 ) jnll_i(i) = -1 * dpois( b_i(i), exp(logchat_i(i)+delta_i(i)), true );
            if( ObsModel_p(p_i(i))==4 ) jnll_i(i) = -1 * dnorm( b_i(i), logchat_i(i), exp(logsigma_pz(p_i(i),0)), true );
        }
    }
        
    // Probability of overdispersion
    for(int i=0; i<n_i; i++){
        if( ObsModel_p(p_i(i))==3 ){
            jnll_i(i) -= dnorm( delta_i(i), Type(0.0), exp(logsigma_pz(p_i(i),0)), true);
        }
    }
        
    // Combine NLL
    jnll_comp(3) = jnll_i.sum();
    jnll = jnll_comp.sum(); 
        
    // Calculate indices
    array<Type> Index_ktp(n_k, n_t, n_p*2);
    matrix<Type> Index_tp(n_t, n_p*2);
    matrix<Type> ln_Index_tp(n_t, n_p*2);
    Index_tp.setZero();
    for(int t=0; t<n_t; t++){
    for(int p=0; p<2*n_p; p++){
        for(int k=0; k<n_k; k++){
            Index_ktp(k,t,p) = d_pkt(p,k,t)  * a_k(k);
            Index_tp(t,p) += Index_ktp(k,t,p);
        }
        ln_Index_tp(t,p) = log( Index_tp(t,p) );
    }}
    REPORT( Index_tp );
    ADREPORT( Index_tp );
    ADREPORT( ln_Index_tp );
        
    // Parameters
    REPORT( R_mu );
    REPORT( phi_p );
    REPORT( logkappa_z );
    REPORT( Hinput_z );
    REPORT( logsigma_pz );    
        
    // Spatial field summaries
    REPORT( Range_pz );
    REPORT( Cov_pp );
    REPORT( logtauE_p );
    REPORT( logtauR );
    REPORT( MargSigmaR );
    REPORT( L_pj );
    // Fields
    REPORT( R_kt );
    REPORT( d_pkt );
    REPORT( dhat1_pkt );
    REPORT( dhat2_pkt );
    REPORT( Rinput_kt );
    // Harvest
    REPORT( F_male_pkt );
    REPORT( chat_pkt );
    REPORT( dpred_pkt );
    // Sanity checks
    REPORT( L_jp );
    REPORT( M_PI );
    REPORT( logchat_i );
    // Derived summaries
    // Objective function components
    REPORT( jnll );
    REPORT( jnll_comp );
    REPORT( jnll_i );
    REPORT( Q_spatiotemporal1 );
    REPORT( Q_spatiotemporal2 );

    return jnll;
 }
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
