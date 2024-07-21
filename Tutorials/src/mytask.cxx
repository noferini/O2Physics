// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Both tasks, ATask and BTask create two histograms. But whereas in
///        the first case (ATask) the histograms are not saved to file, this
///        happens automatically if OutputObj<TH1F> is used to create a
///        histogram. By default the histogram is saved to file
///        AnalysisResults.root. HistogramRegistry is yet an other possibility
///        to deal with histograms. See tutorial example histogramRegistery.cxx
///        for details.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "CommonConstants/LHCConstants.h"
#include "TOFBase/EventTimeMaker.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

namespace o2::aod {

namespace pids {
  DECLARE_SOA_INDEX_COLUMN(Collision, collision);
  DECLARE_SOA_COLUMN(Pt, pt, float);
  DECLARE_SOA_COLUMN(Eta, eta, float);
  DECLARE_SOA_COLUMN(Phi, phi, float);
  DECLARE_SOA_COLUMN(TOFChi2, tofchi2, float);
  DECLARE_SOA_COLUMN(Length, length, float);
  DECLARE_SOA_COLUMN(TOFSignal, tofsignal, float);
  DECLARE_SOA_COLUMN(TOFexpPi, tofexppi, float);
  DECLARE_SOA_COLUMN(TOFexpKa, tofexpka, float);
  DECLARE_SOA_COLUMN(TOFexpPr, tofexppr, float);
  DECLARE_SOA_COLUMN(TOFexpDe, tofexpde, float);
  DECLARE_SOA_COLUMN(TOFexpHe, tofexphe, float);
  DECLARE_SOA_COLUMN(TOFresPi, tofrespi, float);
  DECLARE_SOA_COLUMN(TOFresKa, tofreska, float);
  DECLARE_SOA_COLUMN(TOFresPr, tofrespr, float);
  DECLARE_SOA_COLUMN(TOFresDe, tofresde, float);
  DECLARE_SOA_COLUMN(TOFresHe, tofreshe, float);
  DECLARE_SOA_COLUMN(DeDx, dedx, float);
  DECLARE_SOA_COLUMN(Ptpc, ptpc, float);
  DECLARE_SOA_COLUMN(FT0Qual, ft0qual, int);
}
namespace eventtimes {
  DECLARE_SOA_INDEX_COLUMN(Collision, collision);
  DECLARE_SOA_COLUMN(FT0A, ft0a, float);
  DECLARE_SOA_COLUMN(FT0C, ft0c, float);
  DECLARE_SOA_COLUMN(PVZ, pvz, float);
  DECLARE_SOA_COLUMN(TOFET, tofet, float);
  DECLARE_SOA_COLUMN(TOFETerr, tofeterr, float);
  DECLARE_SOA_COLUMN(TOFETmult, tofetmult, int);
  DECLARE_SOA_COLUMN(TOFETmultOr, tofetmultor, int);
}
}

namespace o2::aod {
DECLARE_SOA_TABLE(MyPID, "AOD", "MYPID",
                  pids::CollisionId, pids::Pt, pids::Eta, pids::Phi, pids::TOFChi2, pids::Length, pids::TOFSignal, pids::TOFexpPi, pids::TOFexpKa, pids::TOFexpPr, pids::TOFexpDe, pids::TOFexpHe, pids::TOFresPi, pids::TOFresKa, pids::TOFresPr, pids::TOFresDe, pids::TOFresHe, pids::DeDx, pids::Ptpc, pids::FT0Qual);
DECLARE_SOA_TABLE(EventTime, "AOD", "EVENTTIME",
                  eventtimes::CollisionId, eventtimes::FT0A, eventtimes::FT0C, eventtimes::PVZ, eventtimes::TOFET, eventtimes::TOFETerr, eventtimes::TOFETmult, eventtimes::TOFETmultOr);
}

using namespace o2;
using namespace o2::framework;

struct Mytask {
  Produces<o2::aod::MyPID> pids;
  Produces<o2::aod::EventTime> eventtimes;

  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<o2::soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal>::iterator, pid>;

  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;

  //  void process(aod::Collision const& collision, o2::soa::Join<aod::Tracks, aod::TracksExtra> const& tracks) {
  void process(o2::soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>::iterator const& collision, o2::soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal> const& tracks, aod::FT0s const& /*ft0s*/) {
    mRespParamsV2.setParameter(4, 68); // TOF resolution
    if(! collision.has_foundFT0()){
      return;
    }
    auto ft0 = collision.foundFT0();
    if(fabs(ft0.timeA()) > 1 || fabs(ft0.timeC()) > 1){
      return;
    }
    float zpv = collision.posZ();
    float t0a = ft0.timeA()*1000 + 32.356410 * zpv; // I need to assume a speed 3% higher than light, why? Is t0a-t0c distance and z-PV set properly?
    float t0c = ft0.timeC()*1000 - 32.356410 * zpv;
    float t0ac = (t0a+t0c)*0.5;
    bool ft0qual = fabs(t0a-t0c) < 100;
    //    printf("FT0A=%f ps,  FT0C=%f ps\n",t0a, t0c);
    //    printf("TA,TC,Z= %f %f %f\n",t0a,t0c,zpv);

    int clId = -1;

    std::vector<tof::eventTimeTrackTest> toftracks;
    constexpr auto responsePi = ResponseImplementation<track::PID::Pion>();
    constexpr auto responseKa = ResponseImplementation<track::PID::Kaon>();
    constexpr auto responsePr = ResponseImplementation<track::PID::Proton>();
    constexpr auto responseDe = ResponseImplementation<track::PID::Deuteron>();
    constexpr auto responseHe = ResponseImplementation<track::PID::Helium3>();

    int ntofinput = 0;
    for (auto& track : tracks) {
      float dedx = track.tpcSignal();
      if(track.tofChi2() < -1 || dedx > 1E9 || dedx < 5 || track.itsChi2NCl() < 0){ // keep only TOF tracks and TPC and ITS signals ok
	continue;
      } else {
	float ptsig = 1./track.signed1Pt();
	clId = track.collisionId();
	tof::eventTimeTrackTest trk;
	
	trk.mTOFChi2 = track.tofChi2();
	trk.mP = fabs(ptsig*1.7);
	trk.expTimes[0] = responsePi.GetCorrectedExpectedSignal(mRespParamsV2, track);
	trk.expTimes[1] = responseKa.GetCorrectedExpectedSignal(mRespParamsV2, track);
	trk.expTimes[2] = responsePr.GetCorrectedExpectedSignal(mRespParamsV2, track);
	float expDe = responseDe.GetCorrectedExpectedSignal(mRespParamsV2, track);
	float expHe = responseHe.GetCorrectedExpectedSignal(mRespParamsV2, track);
	trk.expSigma[0] = responsePi.GetExpectedSigmaTracking(mRespParamsV2, track);
	trk.expSigma[1] = responseKa.GetExpectedSigmaTracking(mRespParamsV2, track);
	trk.expSigma[2] = responsePr.GetExpectedSigmaTracking(mRespParamsV2, track);
	float resDe = responseDe.GetExpectedSigmaTracking(mRespParamsV2, track);
	float resHe = responseHe.GetExpectedSigmaTracking(mRespParamsV2, track);
	trk.mSignal = track.tofSignal();

	if(trk.mP < 2){ // current cut on DummyFilter (to be replaced)
	  ntofinput++;
	}
	
	toftracks.push_back(trk);

	trk.mSignal -= t0ac; // track table is stored after subtracting FT0-AC

        pids(clId, ptsig, track.eta(), track.phi(), track.tofChi2(), track.length(), trk.mSignal, trk.expTimes[0], trk.expTimes[1], trk.expTimes[2], expDe, expHe, trk.expSigma[0], trk.expSigma[1], trk.expSigma[2], resDe, resHe, dedx, track.tpcInnerParam(), ft0qual);
      }
    }
    if(clId > -1){
      auto evtime = tof::evTimeMaker<std::vector<tof::eventTimeTrackTest>, tof::eventTimeTrackTest, tof::filterDummy>(toftracks);
      eventtimes(clId,t0a,t0c,zpv,evtime.mEventTime,evtime.mEventTimeError,evtime.mEventTimeMultiplicity,ntofinput);
    }    
  }
};


struct HistRegistry {

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {{"tofET", "tof Event Time; TOF used tracks; t_{TOF}^{0} - t_{FT0-AC}^{0} (ps)", {HistType::kTH2F, {{40, 0, 40}, {100, -1000.0, 1000.0}}}},
     {"tofETresTh", "tof Event Time Theoretical resolution; TOF used tracks; #sigma (ps)", {HistType::kTProfile, {{40, 0, 40}}}},
     {"effES", "Event section efficiency; N contributors; #varepsilon", {HistType::kTProfile, {{200, 0, 200}}}},
     {"tofPIDsig", "tofPID N#sigma;p_{T} (GeV/#it{c}); N#sigma^{#pi}", {HistType::kTH2F, {{200, 0, 20}, {100, -10.0, 10.0}}}},
     {"tofPID", "tofPID;p_{T} (GeV/#it{c}); t_{TOF} - t_{FT0-AC}^{0} - t_{exp}^{#pi} (ps)", {HistType::kTH2F, {{200, 0, 20}, {100, -1000.0, 1000.0}}}}}};

  
  //  void process(aod::Collision const& collision, aod::MyPID const& tracks)
  void process(aod::Collision const& coll, aod::EventTime const& collisions, aod::MyPID const& tracks)
  {
    int ncontr = coll.numContrib();
    registry.get<TProfile>(HIST("effES"))->Fill(ncontr, collisions.size()>0);
    
    if(!collisions.size()){
      return;
    }
    
    for (auto& collision : collisions) {
      if(fabs(collision.ft0a() - collision.ft0c()) < 100){
	float t0ac = (collision.ft0a() + collision.ft0c())*0.5;
	int mult = collision.tofetmult();
	if(mult > 39){
	  mult = 39;
	}
	registry.get<TH2>(HIST("tofET"))->Fill(mult, collision.tofet() - t0ac);
	float res= collision.tofeterr();
	res = sqrt(res*res + 15*15);
	registry.get<TProfile>(HIST("tofETresTh"))->Fill(mult, res);
      }
    }
    for (auto& track : tracks) {
      if(track.ft0qual()){
	float sigmaPi = track.tofrespi();
	sigmaPi = sqrt(sigmaPi*sigmaPi + 15*15);
	registry.get<TH2>(HIST("tofPIDsig"))->Fill(fabs(track.pt()), (track.tofsignal() - track.tofexppi())/sigmaPi);
	registry.get<TH2>(HIST("tofPID"))->Fill(fabs(track.pt()), track.tofsignal() - track.tofexppi());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Mytask>(cfgc),
    adaptAnalysisTask<HistRegistry>(cfgc),
  };
}
