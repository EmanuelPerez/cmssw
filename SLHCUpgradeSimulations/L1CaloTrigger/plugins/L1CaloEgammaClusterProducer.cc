
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/L1CaloAlgoBase.h"

#include "SimDataFormats/SLHC/interface/L1CaloClusterWithSeed.h"
#include "SimDataFormats/SLHC/interface/L1CaloClusterWithSeedFwd.h"

#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"


class L1CaloEgammaClusterProducer:public L1CaloAlgoBase < l1slhc::L1CaloClusterWithSeedCollection , l1slhc::L1CaloClusterWithSeedCollection  > 
{
  public:
	L1CaloEgammaClusterProducer( const edm::ParameterSet & );
	 ~L1CaloEgammaClusterProducer(  );

//	void initialize(  );

	void algorithm( const int &, const int & );


};

L1CaloEgammaClusterProducer::L1CaloEgammaClusterProducer( const edm::ParameterSet & aConfig ):
L1CaloAlgoBase < l1slhc::L1CaloClusterWithSeedCollection , l1slhc::L1CaloClusterWithSeedCollection > ( aConfig )
{
}

L1CaloEgammaClusterProducer::~L1CaloEgammaClusterProducer(  )
{
}

/*
void L1CaloEgammaClusterProducer::initialize(  )
{
}
*/

void L1CaloEgammaClusterProducer::algorithm( const int &aEta, const int &aPhi )
{

    // Look if there is a cluster here
    l1slhc::L1CaloClusterWithSeedCollection::const_iterator lClusterItr = fetch( aEta, aPhi );
    if ( lClusterItr != mInputCollection->end(  ) )
    {
        // Check if the proto-cluster has been flagged as e/g
        if( lClusterItr->isEGamma() )
        {
            l1slhc::L1CaloClusterWithSeed lEgammaCluster( *lClusterItr );
            int posEtaPlus  = lEgammaCluster.hasConstituent(1, 0);
            int posEtaMinus = lEgammaCluster.hasConstituent(-1, 0);

            int posEtaPlusPhiPlus   = lEgammaCluster.hasConstituent(1, 1);
            int posEtaMinusPhiPlus  = lEgammaCluster.hasConstituent(-1, 1);
            int posEtaMinusPhiMinus = lEgammaCluster.hasConstituent(-1, -1);
            int posEtaPlusPhiMinus  = lEgammaCluster.hasConstituent(1, -1);

            int EetaPlus  = (posEtaPlus!=-1  ? lEgammaCluster.constituentEmEt(1, 0) : 0);
            int EetaMinus = (posEtaMinus!=-1 ? lEgammaCluster.constituentEmEt(-1, 0) : 0);

            int EetaPlusPhiPlus   = (posEtaPlusPhiPlus!=-1    ? lEgammaCluster.constituentEmEt(1, 1) : 0);
            int EetaMinusPhiPlus  = (posEtaMinusPhiPlus!=-1   ? lEgammaCluster.constituentEmEt(-1, 1) : 0);
            int EetaMinusPhiMinus = (posEtaMinusPhiMinus!=-1  ? lEgammaCluster.constituentEmEt(-1, -1) : 0);
            int EetaPlusPhiMinus  = (posEtaPlusPhiMinus!=-1   ? lEgammaCluster.constituentEmEt(1, -1) : 0);

            bool keepEtaPlus = false;
            bool keepEtaMinus = false;

            // First choose to remove or keep towers at eta +/- 1
            if(EetaPlus>EetaMinus) // -> remove eta-1
            {
                keepEtaPlus = true;
            }
            else if(EetaPlus<EetaMinus) // -> remove eta+1
            {
                keepEtaMinus = true;
            }
            else // eta+1 = eta-1 -> look at the corners of the 3x3 proto-cluster
            {
                int EtotEtaPlus  = EetaPlusPhiPlus + EetaPlusPhiMinus;
                int EtotEtaMinus = EetaMinusPhiPlus + EetaMinusPhiMinus;
                if(EtotEtaPlus>EtotEtaMinus) // -> remove eta-1
                {
                    keepEtaPlus = true;
                }
                else if(EtotEtaPlus<EtotEtaMinus) // -> remove eta+1
                {
                    keepEtaMinus = true;
                }
                else // -> keep both eta+1 and eta-1
                {
                    keepEtaPlus = true;
                    keepEtaMinus = true;
                }
            }
            if(!keepEtaPlus)
            {
                lEgammaCluster.removeConstituent(1, -1);
                lEgammaCluster.removeConstituent(1, 0);
                lEgammaCluster.removeConstituent(1, 1);
            }
            if(!keepEtaMinus)
            {
                lEgammaCluster.removeConstituent(-1, -1);
                lEgammaCluster.removeConstituent(-1, 0);
                lEgammaCluster.removeConstituent(-1, 1);
            }
            // check if the towers in the corners are kept or not (they have to be adjacent to a non-zero tower)
            if(keepEtaPlus)
            {
                if(lEgammaCluster.hasConstituent(0, 1)==-1 && lEgammaCluster.hasConstituent(1, 0)==-1)
                {
                    lEgammaCluster.removeConstituent(1, 1);
                }
                if(lEgammaCluster.hasConstituent(0, -1)==-1 && lEgammaCluster.hasConstituent(1, 0)==-1)
                {
                    lEgammaCluster.removeConstituent(1, -1);
                }
            }
            if(keepEtaMinus)
            {
                if(lEgammaCluster.hasConstituent(0, 1)==-1 && lEgammaCluster.hasConstituent(-1, 0)==-1)
                {
                    lEgammaCluster.removeConstituent(-1, 1);
                }
                if(lEgammaCluster.hasConstituent(0, -1)==-1 && lEgammaCluster.hasConstituent(-1, 0)==-1)
                {
                    lEgammaCluster.removeConstituent(-1, -1);
                }
            }


            int lIndex = mCaloTriggerSetup->getBin( aEta, aPhi );
            std::pair < int, int >lEtaPhi = mCaloTriggerSetup->getTowerEtaPhi( lIndex );
            mOutputCollection->insert( lEtaPhi.first , lEtaPhi.second , lEgammaCluster );
        }
    }
}



DEFINE_EDM_PLUGIN (edm::MakerPluginFactory,edm::WorkerMaker<L1CaloEgammaClusterProducer>,"L1CaloEgammaClusterProducer");
DEFINE_FWK_PSET_DESC_FILLER(L1CaloEgammaClusterProducer);

