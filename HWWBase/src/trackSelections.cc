#include "Math/VectorUtil.h"
#include <math.h>
#include "HWWValidation/HWWBase/interface/HWW.h"
#include "HWWValidation/HWWBase/interface/trackSelections.h"

// return a pair of d0, d0err of a ctf track with respect to a primary vertex
std::pair<double, double> trks_d0_pv (int itrk, int ipv)
{
    // assume the layout of the covariance matrix is (Vxx, Vxy, Vxz)
    //						      (Vyx, Vyy, ...)
    const double bx  = HWWVal::vtxs_position().at(ipv).x()   ;
    const double by  = HWWVal::vtxs_position().at(ipv).y()   ;
    const double vxx = HWWVal::vtxs_covMatrix().at(ipv).at(0);
    const double vxy = HWWVal::vtxs_covMatrix().at(ipv).at(1);
    const double vyy = HWWVal::vtxs_covMatrix().at(ipv).at(4);

    const double phi      = HWWVal::trks_trk_p4().at(itrk).phi();
    const double d0vtx    = HWWVal::trks_d0().at(itrk) - bx * sin(phi) + by * cos(phi);
    const double d0err    = HWWVal::trks_d0Err().at(itrk);
    const double phierr   = HWWVal::trks_phiErr().at(itrk);
    const double d0phicov = HWWVal::trks_d0phiCov().at(itrk);

    // we will let the optimizer take care of subexpression
    // elimination for this one...
    const double d0err2vtx = d0err * d0err 
        - 2 * (bx * cos(phi) + by * sin(phi)) * d0phicov
        + (bx * cos(phi) + by * sin(phi)) * (bx * cos(phi) + by * sin(phi)) * phierr * phierr
        + sin(phi) * sin(phi) * vxx + cos(phi) * cos(phi) * vyy
        - 2 * sin(phi) * cos(phi) * vxy;
    if (d0err2vtx >= 0) 
        return std::pair<double, double>(d0vtx, sqrt(d0err2vtx));

    std::cerr << "Oh no!  sigma^2(d0corr) < 0!" << std::endl;
    return std::pair<double, double>(d0vtx, -sqrt(-d0err2vtx));
}

// return a pair of d0, d0err of a gsf track with respect to a primary vertex
std::pair<double , double> gsftrks_d0_pv (int itrk, int ipv)
{
    // assume the layout of the covariance matrix is (Vxx, Vxy, Vxz)
    //						      (Vyx, Vyy, ...)
    const double bx  = HWWVal::vtxs_position().at(ipv).x()   ;
    const double by  = HWWVal::vtxs_position().at(ipv).y()   ;
    const double vxx = HWWVal::vtxs_covMatrix().at(ipv).at(0);
    const double vxy = HWWVal::vtxs_covMatrix().at(ipv).at(1);
    const double vyy = HWWVal::vtxs_covMatrix().at(ipv).at(4);

    const double phi      = HWWVal::gsftrks_p4().at(itrk).phi();
    const double d0vtx    = HWWVal::gsftrks_d0().at(itrk) - bx * sin(phi) + by * cos(phi);
    const double d0err    = HWWVal::gsftrks_d0Err().at(itrk);
    const double phierr   = HWWVal::gsftrks_phiErr().at(itrk);
    const double d0phicov = HWWVal::gsftrks_d0phiCov().at(itrk);

    // we will let the optimizer take care of subexpression
    // elimination for this one...
    const double d0err2vtx = d0err * d0err 
        - 2 * (bx * cos(phi) + by * sin(phi)) * d0phicov
        + (bx * cos(phi) + by * sin(phi)) * (bx * cos(phi) + by * sin(phi)) * phierr * phierr
        + sin(phi) * sin(phi) * vxx + cos(phi) * cos(phi) * vyy
        - 2 * sin(phi) * cos(phi) * vxy;
    if (d0err2vtx >= 0) 
        return std::pair<double, double>(d0vtx, sqrt(d0err2vtx));

    std::cerr << "Oh no!  sigma^2(d0corr) < 0!" << std::endl;
    return std::pair<double, double>(d0vtx, -sqrt(-d0err2vtx));
}

// return a pair of dz, dzerr of a ctf track with respect to a primary vertex
std::pair<double, double> trks_dz_pv (int itrk, int ipv)
{


    LorentzVector pv = HWWVal::vtxs_position().at(ipv);
    double pvxErr    = HWWVal::vtxs_xError().at(ipv)  ;
    double pvyErr    = HWWVal::vtxs_yError().at(ipv)  ;
    double pvzErr    = HWWVal::vtxs_zError().at(ipv)  ;
    double phi        = HWWVal::trks_trk_p4().at(itrk).phi();
    double theta      = HWWVal::trks_trk_p4().at(itrk).theta();
    double ddzdpvx    = cos(phi)*1./tan(theta);
    double ddzdpvy    = sin(phi)*1./tan(theta);
    double ddzdphi    = -1*pv.x()*sin(phi)*1./tan(theta) + pv.y()*cos(phi)*1./tan(theta);
    double ddzdtheta  = -1*1/sin(theta)*1/sin(theta) * (pv.x()*cos(phi) + pv.y()*sin(phi));

    ddzdpvx   *= ddzdpvx;
    ddzdpvy   *= ddzdpvy;
    ddzdphi   *= ddzdphi;
    ddzdtheta *= ddzdtheta;

    double z0Err    = HWWVal::trks_z0Err().at(itrk);
    double phiErr   = HWWVal::trks_phiErr().at(itrk);
    double thetaErr = HWWVal::trks_etaErr().at(itrk)*sin(theta);

    z0Err    *= z0Err;
    phiErr   *= phiErr;
    thetaErr *= thetaErr;
    pvxErr   *= pvxErr;
    pvyErr   *= pvyErr;
    pvzErr   *= pvzErr;

    double value = HWWVal::trks_z0().at(itrk) - pv.z() + (pv.x()*cos(phi) + pv.y()*sin(phi) )*1./tan(theta);

    //note that the error does not account for correlations since we do not store the track covariance matrix
    double error = sqrt(z0Err + pvzErr + ddzdpvx*pvxErr + ddzdpvy*pvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

    return std::pair<double, double>(value, error);
}

std::pair<double, double> gsftrks_dz_pv (int itrk, int ipv)
{
    LorentzVector pv = HWWVal::vtxs_position().at(ipv);
    double pvxErr    = HWWVal::vtxs_xError().at(ipv)  ;
    double pvyErr    = HWWVal::vtxs_yError().at(ipv)  ;
    double pvzErr    = HWWVal::vtxs_zError().at(ipv)  ;

    //LorentzVector tkp = HWWVal::gsftrks_p4().at(itrk);
    //LorentzVector tkv = HWWVal::gsftrks_vertex_p4().at(itrk);

    double phi   = HWWVal::gsftrks_p4().at(itrk).phi();
    double theta = HWWVal::gsftrks_p4().at(itrk).theta();

    double ddzdpvx   = cos(phi)*1./tan(theta);
    double ddzdpvy   = sin(phi)*1./tan(theta);
    double ddzdphi   = -1*pv.x()*sin(phi)*1./tan(theta) + pv.y()*cos(phi)*1./tan(theta);
    double ddzdtheta = -1*1/sin(theta)*1/sin(theta) * (pv.x()*cos(phi) + pv.y()*sin(phi));

    ddzdpvx   *= ddzdpvx;
    ddzdpvy   *= ddzdpvy;
    ddzdphi   *= ddzdphi;
    ddzdtheta *= ddzdtheta;

    double z0Err    = HWWVal::gsftrks_z0Err().at(itrk);
    double phiErr   = HWWVal::gsftrks_phiErr().at(itrk);
    double thetaErr = HWWVal::gsftrks_etaErr().at(itrk)*sin(theta);

    z0Err    *= z0Err;
    phiErr   *= phiErr;
    thetaErr *= thetaErr;
    pvxErr   *= pvxErr;
    pvyErr   *= pvyErr;
    pvzErr   *= pvzErr;

    double value = HWWVal::gsftrks_z0().at(itrk) - pv.z() + (pv.x()*cos(phi) + pv.y()*sin(phi) )*1./tan(theta);

    //note that the error does not account for correlations since we do not store the track covariance matrix
    double error = sqrt(z0Err + pvzErr + ddzdpvx*pvxErr + ddzdpvy*pvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

    return std::pair<double, double>(value, error);
}
