//
//  ThreeBodyUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/10/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef ThreeBodyUnitTest_h
#define ThreeBodyUnitTest_h

//#include "ThreeBody.h"

class TBDUTest : public ::testing::Test
{
  public :
  
  protected :
  
  int vals_;
  vector< Molecule > mol_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -8.0 ), Pt( 0, 0, 0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -7.0 ), Pt( 0, 0, 1 ) };
    double cg[2] = { -5.0, -5.0};
    double rd[2] = {  3.6,  3.6};
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 3; vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0]=cg[molInd]; vdW[0]=0.0; posCharges[0]=cgPos[molInd];
      charges[1]=cg[molInd]; vdW[1]=0.0; posCharges[1]=pos[molInd]+Pt(1,0,0);
      charges[2]=cg[molInd]; vdW[2]=0.0; posCharges[2]=pos[molInd]+Pt(0,1,0);
      
      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd]);
      mol_.push_back( molNew );
    }
  } // end SetUp
  
  virtual void TearDown() {}
  
  int di_cut[40] = {0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,1,5,1,6,1,7,1,8,2,3,2,4,
    2,7,2,8,3,7,4,8,5,7,6,8};
  int tr_cut[180] = {0,1,2,0,1,3,0,1,4,0,1,5,0,1,6,0,1,7,0,1,8,0,2,3,0,2,4,0,
    2,5,0,2,6,0,2,7,0,2,8,0,3,4,0,3,5,0,3,6,0,3,7,0,3,8,0,4,5,0,4,6,0,4,7,0,
    4,8,0,5,6,0,5,7,0,5,8,0,6,7,0,6,8,0,7,8,1,2,4,1,2,5,1,2,6,1,2,7,1,2,8,1,
    3,6,1,3,7,1,3,8,1,4,6,1,4,7,1,4,8,1,5,6,1,5,7,1,5,8,1,6,7,1,6,8,1,7,8,2,
    3,4,2,3,7,2,3,8,2,4,7,2,4,8,2,5,7,2,5,8,2,6,8,2,7,8,3,4,8,3,5,7,3,6,8,4,
    5,8,4,6,8,5,6,8};
  
  double nrgDi7[42] = {0.00374634957,0.00374634957,0.00374634957,0.00374634957,0.00447973634,0.00447973634,0.0027528143,0.0027528143,0.0027528143,0.0027528143,0.00447973634,0.00447973634,0.00172273946,0.00172273946,0.00134728982,0.00134728982,0.00115374419,0.00115374419,0.00444541177,0.00444541177,0.00444541177,0.00444541177,0.00444541177,0.00444541177,0.00444541177,0.00444541177,0.00115374419,0.00115374419,0.00134728982,0.00134728982,0.00217378518,0.00217378518,0.00172273946,0.00172273946,0.00101733621,0.00101733621,0.000838510727,0.000838510727,0.00172273946,0.00172273946,0.00217378518,0.00217378518};
  double foxDi7[42] = {-0.000807962967,0.000807949103,0.000808043166,-0.00080805703,0.00103057776,-0.00103057776,0.000492850966,-0.000492850966,-0.000492958814,0.000492958814,-0.0010301614,0.0010301614,0.000305543083,-0.000305543083,0.000212315627,-0.000212314938,0.000162716843,-0.000162716872,-0.000159979555,0.000159990629,0.000160351505,-0.000160340431,-0.000159979555,0.000159990629,0.000160351505,-0.000160340431,-0.000162717523,0.000162717495,-0.000212314255,0.000212314944,5.72541019e-05,-5.72541019e-05,-0.000305550753,0.000305550753,-0.000130726598,0.000130726598,-9.95612241e-05,9.95612241e-05,-0.000305550753,0.000305550753,5.72541019e-05,-5.72541019e-05};
  double foyDi7[42] = {0.000159998372,-0.000160007953,-0.000160332688,0.000160323107,0.000741190279,-0.000741190279,-0.000203759793,0.000203759793,0.000204093191,-0.000204093191,-0.000741592082,0.000741592082,-5.72533021e-05,5.72533021e-05,7.71170715e-05,-7.71188487e-05,-8.13705302e-05,8.13705405e-05,0.00089805392,-0.000898037435,-0.000898131789,0.000898148273,0.00089805392,-0.000898037435,-0.000898131789,0.000898148273,8.13703902e-05,-8.13703798e-05,-7.71206381e-05,7.71188609e-05,-0.000456740131,0.000456740131,5.72435547e-05,-5.72435547e-05,-0.000116931109,0.000116931109,8.57654403e-05,-8.57654403e-05,5.72435547e-05,-5.72435547e-05,-0.000456740131,0.000456740131};
  double fozDi7[42] = {-0.000506776518,0.000506753299,-0.000506418574,0.000506441793,7.16377524e-20,-7.46266015e-20,-1.43021039e-21,-8.88038268e-22,-2.21133667e-20,2.17508409e-20,7.49701126e-20,-7.4000701e-20,1.5880295e-20,-1.59270034e-20,9.38704413e-05,-9.38708364e-05,5.87462492e-05,-5.87462663e-05,0.00088582343,-0.000885791755,0.000885496733,-0.000885528408,0.00088582343,-0.000885791755,0.000885496733,-0.000885528408,5.87473919e-05,-5.87473748e-05,9.38712374e-05,-9.38708423e-05,3.09823316e-20,-3.03352355e-20,1.0158897e-20,-9.74863537e-21,1.066771e-20,-1.09165638e-20,5.84089332e-21,-5.59879369e-21,1.0158897e-20,-9.74863537e-21,3.09823316e-20,-3.03352355e-20};
  double toxDi7[42] = {-0.000906150555,0.000377358096,-0.000377404696,0.000905485521,4.93630391e-20,-8.1526567e-20,2.62174978e-21,-4.79160869e-20,7.4319912e-20,2.03877719e-20,1.11238565e-19,-3.48104141e-20,1.22018548e-20,-2.99124227e-20,6.89758094e-05,-0.000134069782,4.82884177e-05,-0.000119488089,0.000634160571,-0.000757382277,0.000757801057,-0.000633518371,0.000634160571,-0.000757382277,0.000757801057,-0.000633518371,0.000119490758,-4.82878711e-05,0.000134072307,-6.89746482e-05,3.15849649e-20,-2.05587934e-20,3.82836688e-20,-7.24124716e-21,1.44088116e-20,-7.1159518e-21,1.23169905e-20,-4.16937876e-21,3.82836688e-20,-7.24124716e-21,3.15849649e-20,-2.05587934e-20};
  double toyDi7[42] = {-0.00133534311,0.000177174686,-0.000177182803,0.00133436966,4.06717561e-20,-6.84131164e-20,2.08481361e-20,-1.75731112e-19,2.06277148e-19,-1.00306794e-21,9.78659106e-20,-2.37737175e-20,6.83723133e-21,-4.94224894e-20,5.09938433e-05,-0.00015139619,3.19276614e-05,-0.000205738376,0.000896530657,-0.000367529568,0.000367784798,-0.000895656993,0.000896530657,-0.000367529568,0.000367784798,-0.000895656993,0.000205742361,-3.19272982e-05,0.000151399075,-5.09928417e-05,1.67860755e-20,-3.41623031e-20,6.14391645e-20,-4.36962312e-21,1.37104825e-20,-6.16598767e-21,2.14278074e-20,-2.71522138e-21,6.14391645e-20,-4.36962312e-21,1.67860755e-20,-3.41623031e-20};
  double tozDi7[42] = {0.00105840172,-0.000515052782,-0.000515281061,0.00105758758,-0.00081135645,0.00209180434,-0.000225886627,-0.00105546577,-0.00105595353,-0.000225632453,0.00209319661,-0.000810706112,-0.000218108621,0.000612126699,-0.000196752,0.000428140231,-9.10073433e-05,6.04163365e-05,-0.000758654292,0.000215427512,0.000215533629,-0.000757962307,-0.000758654292,0.000215427512,0.000215533629,-0.000757962307,6.04183418e-05,-9.10059264e-05,0.000428148367,-0.000196748441,0.000173473207,-0.000567487194,0.000612173727,-0.000218101713,0.0002966643,-0.000171787902,-8.70114087e-05,-3.78632749e-05,0.000612173727,-0.000218101713,0.000173473207,-0.000567487194};
  
  double nrgTr7[105] = {0.00752852317,0.00552188109,0.00554390432,0.00816977454,0.00508988575,0.00583115888,0.00652335486,0.00494247114,0.00397574427,0.00664801467,0.00830524722,0.00730739864,0.00838994308,0.00833843474,0.00908194999,0.00839053831,0.00837226231,0.00908315461,0.00661765227,0.00832110344,0.00734289152,0.00644340722,0.00491300148,0.00389938019,0.00825425964,0.00515119149,0.00587830385,0.00733790588,0.00675589558,0.00505549249,0.00718926084,0.00622553924,0.00447929701,0.00889859443,0.00549375486,0.00549447492,0.00547459691,0.00361994036,0.00359668634,0.00729883155,0.00458503787,0.00628267683,0.00724405309,0.0049548023,0.00669828041,0.00307946481,0.00616818279,0.00580708245,0.00289735628,0.00618466692,0.0056341578,0.00619135304,0.00290896144,0.00562723642,0.00620526767,0.00310971191,0.00583561173,0.00251377808,0.00353470625,0.00334357917,0.00583561173,0.00310971191,0.00620526767,0.00581402791,0.00238479123,0.00548857291,0.00561463962,0.00201055834,0.00530143972,0.00562723642,0.00290896144,0.00619135304,0.00896235128,0.00669965822,0.00672347561,0.00896235128,0.00669965822,0.00672347561,0.0056341578,0.00618466692,0.00289735628,0.00580193423,0.00547666531,0.00237578711,0.00560635808,0.00530147049,0.00200208247,0.00580708245,0.00616818279,0.00307946481,0.00251751876,0.00334155495,0.00353637908,0.0039072681,0.00302989766,0.00257318915,0.00319765031,0.00389819601,0.00274595524,0.00275972956,0.00391385375,0.00321110396,0.00257731654,0.00302367747,0.00390582523};
  double foxTr7[105] = {-7.04344471e-06,0.00113645545,-0.00113954955,0.000204226836,0.00102748947,-0.00124222454,-0.000314055416,0.00098404845,-0.00068391081,-0.00132353319,0.000645670654,0.000672401747,-0.00188699856,0.00100506424,0.000881047068,0.00190776416,-0.00099098685,-0.00090340541,0.00135505886,-0.000680414825,-0.00066059424,0.000314994665,-0.000977621353,0.000654275705,-0.000230081995,-0.00103935294,0.00126607754,0.00156505156,-0.00099259825,-0.000577718237,0.000528754372,-0.00134714221,0.000809030553,-2.46931289e-05,-0.00115474057,0.00117087283,8.91545198e-06,-0.000608989298,0.000591028751,-0.000541032045,-0.000837205566,0.00136353907,-0.00147750399,0.000542959032,0.000971844087,0.000520036108,-0.000449450437,-5.89964665e-05,0.000474179409,-0.000142878668,-0.0003244911,0.000141512576,-0.000478602631,0.000329510092,0.000471606365,-0.000530208934,5.51963937e-05,0.000378404149,-0.000157136081,-0.000220416776,5.55715234e-05,-0.000530213278,0.000471265375,0.000374180746,-0.000347429868,-2.83221907e-05,6.65761109e-07,-0.00026666276,0.000262986043,0.000329873779,-0.000478607514,0.000141244787,-9.55008031e-06,0.000223337048,-0.000214050167,-9.55008031e-06,0.000223337048,-0.000214050167,-0.000324127003,-0.000143346432,0.000474185412,-0.000368387863,2.66830056e-05,0.000345943959,-1.37736816e-06,-0.0002588146,0.000264533853,-5.86747808e-05,-0.000449876607,0.000520043576,-0.000379092154,0.000222246844,0.000155702166,-0.000247905246,-0.000158906291,0.000407986543,-7.62966861e-05,-0.000357078113,0.000437127764,-0.000441795791,0.00036557413,7.5438982e-05,-0.000409604029,0.000159089619,0.000247087232};
  double foyTr7[105] = {8.26829968e-06,-0.000223133873,0.000214250205,0.000890858109,-8.22009237e-05,-0.000824369904,-4.28248441e-05,-0.000248475424,0.000306438924,0.000378378052,0.000762687529,-0.00114035065,-0.000592364382,-0.00109094926,0.00166850902,0.000613511619,0.00107971801,-0.00169178577,-0.000375699717,-0.000738008057,0.00110371063,2.52311238e-05,0.000240260199,-0.000277677478,-0.000886251443,7.27047833e-05,0.00082443825,0.000543537689,-0.00123258376,0.000698999118,0.000873447682,-0.000676887524,-0.000246595308,2.48609283e-05,-0.000868305245,0.000853309033,-1.83400269e-05,0.000308753838,-0.000287257617,-0.000958149226,0.000281435614,0.000694886168,-0.000526614172,-0.000678474195,0.00120935404,2.24095234e-05,0.000948350655,-0.00097727379,-0.000141662543,-0.000836449533,0.000984461045,0.000843846034,0.00014077679,-0.000986229298,-0.000961717543,-2.19020656e-05,0.00098379673,-3.6915241e-06,-0.000536240629,0.000539978995,0.000983640652,-2.19211638e-05,-0.000961714807,-0.000822967313,-0.000197219807,0.00102045312,0.000819205585,0.000170403568,-0.000989068602,-0.000986268852,0.00014077096,0.000843986512,7.81595577e-06,-0.00139031384,0.00139184087,7.81595577e-06,-0.00139031384,0.00139184087,0.000984303589,-0.000836477039,-0.000141647186,0.000822779417,-0.00101956368,0.000196402216,-0.000814609819,0.000984747939,-0.000168919194,-0.000977295179,0.000948431375,2.2419165e-05,4.67043477e-06,-0.000541349177,0.000536314187,-0.000398743571,0.000547615977,-0.000145331532,-0.000575035309,0.000508975025,6.22649346e-05,-6.14793037e-05,-0.000517123222,0.000579165512,0.000144632794,-0.000545248528,0.00039974832};
  double fozTr7[105] = {-0.00100061453,0.000509139055,0.000516042742,-0.000464894483,0.000590858182,-0.000102883114,-0.000500441403,0.000569544612,-6.36983337e-05,-0.00050429462,0.00140892878,-0.000915319701,-0.000534399042,0.00144232419,-0.000908233809,-0.000523039399,0.00143034403,-0.000903479036,-0.000540758422,0.00142567692,-0.000883336082,-0.000488042252,0.000573281927,-6.42032811e-05,-0.000508090678,0.000617753529,-9.94286849e-05,7.21464933e-20,-4.62372138e-20,-3.37363083e-20,4.64072763e-20,-6.49444772e-20,8.62903786e-21,1.43101483e-19,-6.31062612e-20,-8.40089637e-20,-2.27479715e-20,5.8735856e-21,1.59738706e-20,8.19244762e-20,8.21359055e-21,-8.62128284e-20,5.48185203e-20,4.7738157e-20,-1.04116492e-19,9.84592132e-05,0.000873350253,-0.000979775032,6.04096263e-05,0.000884227868,-0.000951801589,0.000886167381,6.16981372e-05,-0.000949687755,0.000892350633,9.57030469e-05,-0.000989455383,0.000154083285,-9.49836558e-05,-5.88939761e-05,0.000989716806,-9.57115228e-05,-0.000892767541,0.000984109418,-9.53552584e-05,-0.000890915238,0.000945252076,-5.94469696e-05,-0.000888093546,0.000949295048,-6.17044867e-05,-0.000885917739,0.00176996022,-0.000893726101,-0.00090222344,0.00176996022,-0.000893726101,-0.00090222344,0.000952065104,-0.000884630286,-6.04109917e-05,0.000974230867,-0.000885826701,-9.29696782e-05,0.000941776705,-0.000887892108,-5.9066681e-05,0.000979343856,-0.000873061987,-9.84509968e-05,0.000154586634,-5.96762395e-05,-9.41504103e-05,4.02606676e-20,-2.49136925e-20,-1.55561478e-20,4.18394871e-20,-2.02576798e-20,-2.04031993e-20,2.12105336e-20,2.1241091e-20,-4.16886019e-20,1.62487619e-20,2.51750442e-20,-3.99787853e-20};
  double toxTr7[105] = {-0.00129559535,0.000375975398,0.000901574273,-0.000909263871,0.000437024193,-0.000144325452,-0.00091073241,0.000426776728,-0.000124569624,-0.000913681474,0.00102640876,-0.000772606791,-0.000932808504,0.001157439,-0.000641483691,-0.000379873695,0.00156113901,-0.000772620472,-0.000396980985,0.0016942048,-0.000638379566,-0.000381373884,0.00101048964,-5.21141581e-05,-0.000383100872,0.00103897081,-7.12086651e-05,5.24860145e-20,-5.21952515e-20,-6.85102776e-20,1.22051247e-19,-4.15458434e-20,1.12685919e-20,1.57376767e-19,-6.58521382e-20,-4.06950049e-20,7.52032274e-20,-3.34977089e-20,1.66681639e-20,1.19505761e-19,-7.45286296e-21,-3.83447363e-20,1.85439566e-19,5.20853397e-20,-4.76508856e-20,7.17212703e-05,0.000635889515,-0.000891261372,4.91057457e-05,0.000763287097,-0.000754253707,0.000635365014,0.000122881521,-0.000806256015,0.000760714903,0.0001358668,-0.000704088138,0.000117929054,-0.00013501996,-0.000119893877,0.000704739418,-0.000135876146,-0.00076034333,0.000829460769,-0.000135057112,-0.000633470813,0.000683703697,-0.000118938285,-0.000757108134,0.000806666729,-0.000122891196,-0.000634741861,0.00140947647,-0.000757972513,-0.000629052276,0.00140947647,-0.000757972513,-0.000629052276,0.000754913625,-0.000762966655,-4.91029272e-05,0.000769977131,-0.000756659366,-6.80201642e-05,0.000881344886,-0.000632987061,-4.85416536e-05,0.000891705177,-0.000635278708,-7.1711394e-05,0.000254642438,-4.89818885e-05,-6.92807533e-05,6.92936508e-20,-8.21713226e-21,-1.15868538e-20,4.53438007e-20,1.77907471e-20,-1.40135001e-20,5.29798364e-20,2.49296147e-20,-2.76034097e-20,5.0536698e-20,2.68985238e-20,-2.79264597e-20};
  double toyTr7[105] = {-0.00152315762,0.000175328226,0.00133224441,-0.00134115446,0.00022298995,-0.00016213384,-0.00133975267,0.000209069927,-0.000217962395,-0.00136658568,0.00108780531,-0.000370284426,-0.0013670962,0.00055261595,-0.000903637177,-0.000179461574,0.00227082867,-0.000377640097,-0.000182854372,0.00173362676,-0.000913025328,-0.000180249814,0.00152054351,-3.35061581e-05,-0.000180442334,0.00147926783,-5.29565499e-05,6.17490787e-20,-5.10519205e-20,-2.11782558e-19,2.4189732e-19,-4.57060305e-21,-3.19777065e-21,1.33683639e-19,-5.36339072e-20,-2.89953311e-20,2.23634273e-19,-1.52145135e-19,-2.84601092e-21,1.1698218e-19,-1.12747709e-19,-2.38231309e-20,2.99180232e-19,2.2342666e-20,-5.48287021e-20,5.26097516e-05,0.000900707588,-0.000521210044,3.23385902e-05,0.000371991737,-0.00110552149,0.000897562011,0.000211794822,-0.00039891405,0.000368411741,0.000154641082,-0.00094755257,8.33362692e-05,-0.000152062557,-0.000206502665,0.000948429019,-0.000154654984,-0.000368192963,0.000419896934,-0.000152128022,-0.000895037598,0.000929487897,-0.000205122987,-0.000366741992,0.000399157903,-0.000211809118,-0.000896705337,0.00127759747,-0.000365472543,-0.000891060177,0.00127759747,-0.000365472543,-0.000891060177,0.00110642943,-0.000371823019,-3.23364305e-05,0.00105098579,-0.000367464186,-5.03126569e-05,0.000577148314,-0.000896347809,-3.20175016e-05,0.000521489532,-0.000899877688,-5.26024715e-05,0.00035892061,-3.22871351e-05,-5.13482546e-05,7.86546427e-20,-1.23718352e-20,-7.19771186e-21,3.07220972e-20,2.73431036e-20,-1.06287116e-20,7.53644789e-20,1.2509905e-20,-4.02501086e-20,8.29027719e-20,1.40077818e-20,-3.86131244e-20};
  double tozTr7[105] = {0.000542865956,-0.000740219467,0.00168725779,0.000245328131,-0.000711135396,0.00249561243,0.000830435155,-0.000609290505,-0.00098219031,-7.86996361e-06,-0.00128428025,-1.80263185e-05,0.00321383108,-0.000318356555,-0.00158653453,-0.00134880702,0.000291861211,0.00235418393,-0.000762940853,0.0013099732,-0.00184380826,-0.00157792327,0.00109909832,-0.000313748317,0.00158460276,0.00150192473,-0.00100948686,-0.00105420688,0.00228840245,-0.00164567452,-0.00187880671,0.00270058779,-0.000448435136,0.00128869882,0.00235756877,-0.00097510346,-0.00128879052,-0.00112202533,-0.000260420416,0.00187216494,-0.00043184766,-0.00103792959,0.00106052237,-4.20846549e-05,-0.00138280479,-0.000415587577,-0.000147729343,0.000653966483,-0.000311507802,0.000837644002,-0.000697755089,-0.00098017317,0.000670512366,0.000119833855,-6.43732921e-06,0.00105161992,-0.000957368758,-0.000289385215,0.000605212159,-0.000508486021,-0.000958068859,0.00105167305,-6.46130996e-06,1.63282045e-05,0.000730622663,-0.00093169982,-0.00085089744,-2.89464284e-05,0.000175661089,0.000119914468,0.00067054646,-0.000979529467,-0.00054283642,0.000393938661,-0.00134099428,-0.00054283642,0.000393938661,-0.00134099428,-0.000698396601,0.000837612854,-0.000311486185,-0.000330280544,0.000517702865,-0.000370124103,0.000280295419,-0.000846697549,-0.00012966784,0.000654241748,-0.000147062074,-0.000415556042,0.000490412704,8.14654394e-05,-0.000765867994,0.00078910536,-0.000655454558,-0.000256756877,0.000473760656,4.40969332e-05,-0.000390540179,0.000913725275,-4.55699498e-05,-0.000741060483,0.000522209447,0.000134254102,-0.00078755881};
  
  double nrg7[7] = {0.0246488808,0.0298503645,0.029007848,0.0210495416,0.0199046973,0.0193436631,0.0222221363};
  double fox7[7] = {-0.000365442356,0.00203345505,-0.00179925288,-0.000930109691,-0.00207224144,0.00232360839,0.000827818726};
  double foy7[7] = {3.27413598e-05,-0.00136602747,0.000434272626,-0.00699222033,0.007163832,-0.00697929868,0.00791708577};
  double foz7[7] = {8.56605538e-21,8.21793859e-19,6.79655944e-19,-3.04997911e-19,-5.13117743e-19,-2.21292688e-19,-6.3579692e-19};
  double tox7[7] = {4.75763543e-19,6.75132516e-19,5.23878167e-19,-2.28292129e-19,6.02574374e-20,-1.89658276e-19,-2.47610253e-20};
  double toy7[7] = {9.44839918e-19,4.75639086e-19,5.36545046e-20,-1.39845808e-20,-3.98370472e-19,2.02269035e-21,-2.90704008e-19};
  double toz7[7] = {0.00305549561,-0.00482526217,0.00235112497,0.00421612423,-0.00417180644,-0.000343075508,-0.00526962743};
  
};

TEST_F(TBDUTest, computeGroups)
{
  vector<Molecule> mol_sing_;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < 9; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_sing_.push_back( molNew );
  }
  
  const int vals = 18;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_sing_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();
  
  EXPECT_TRUE( dim.size() == 36 );
  EXPECT_TRUE( tri.size() == 84 );
  
  EXPECT_TRUE( dim[0][0] == 0 );
  EXPECT_TRUE( dim[0][1] == 1 );
  
  EXPECT_TRUE( dim[35][0] == 7 );
  EXPECT_TRUE( dim[35][1] == 8 );
  
  EXPECT_TRUE( tri[0][0] == 0 );
  EXPECT_TRUE( tri[0][1] == 1 );
  EXPECT_TRUE( tri[0][2] == 2 );
  
  EXPECT_TRUE( tri[83][0] == 6 );
  EXPECT_TRUE( tri[83][1] == 7 );
  EXPECT_TRUE( tri[83][2] == 8 );
}

TEST_F(TBDUTest, computeGroupsCutoff)
{
  vector<Molecule> mol;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < 9; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 18;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest, INTERNAL, 7.5 );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();

  int j;
  for ( j = 0; j < dim.size(); j++)
  {
    EXPECT_TRUE( dim[j][0] == di_cut[j*2] );
    EXPECT_TRUE( dim[j][1] == di_cut[j*2+1] );
  }

  for ( j = 0; j < tri.size(); j++)
  {
    EXPECT_TRUE( tri[j][0] == tr_cut[j*3] );
    EXPECT_TRUE( tri[j][1] == tr_cut[j*3+1] );
    EXPECT_TRUE( tri[j][2] == tr_cut[j*3+2] );
  }
}

TEST_F(TBDUTest, twoBD)
{
  vector<Molecule> mol;
  const int num = 7;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 5.0 ),Pt( -5.0, 0.0, 5.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < num; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 10;
  Constants const_( INTERNAL );
  const_.set_salt_concentration(0.1);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  threeBodTest.solveNmer(2);
  threeBodTest.calcTwoBDEnForTor();
  
  vector<vector<double > > en = threeBodTest.getDiEn();
  vector<vector<Pt > > fr = threeBodTest.getDiFo();
  vector<vector<Pt > > to = threeBodTest.getDiTo();
  
  int j;
  for ( j = 0; j < en.size(); j++)
  {
    if ( en[j][0] != 0) EXPECT_NEAR( en[j][0]/nrgDi7[j*2], 1.0, preclim);
    if ( en[j][1] != 0) EXPECT_NEAR( en[j][1]/nrgDi7[j*2+1], 1.0, preclim);
    
    if ( fr[j][0].x() != 0) EXPECT_NEAR(fr[j][0].x()/foxDi7[j*2],1.0,preclim);
    if ( fr[j][1].x() != 0) EXPECT_NEAR(fr[j][1].x()/foxDi7[j*2+1],1.0,preclim);
    
    if ( fr[j][0].y() != 0) EXPECT_NEAR(fr[j][0].y()/foyDi7[j*2],1.0,preclim);
    if ( fr[j][1].y() != 0) EXPECT_NEAR(fr[j][1].y()/foyDi7[j*2+1],1.0,preclim);
    
    if ( fr[j][0].z() != 0) EXPECT_NEAR(fr[j][0].z()/fozDi7[j*2],1.0,preclim);
    if ( fr[j][1].z() != 0) EXPECT_NEAR(fr[j][1].z()/fozDi7[j*2+1],1.0,preclim);
    
    if ( to[j][0].x() != 0) EXPECT_NEAR(to[j][0].x()/toxDi7[j*2],1.0,preclim);
    if ( to[j][1].x() != 0) EXPECT_NEAR(to[j][1].x()/toxDi7[j*2+1],1.0,preclim);
    
    if ( to[j][0].y() != 0) EXPECT_NEAR(to[j][0].y()/toyDi7[j*2],1.0,preclim);
    if ( to[j][1].y() != 0) EXPECT_NEAR(to[j][1].y()/toyDi7[j*2+1],1.0,preclim);
    
    if ( to[j][0].z() != 0) EXPECT_NEAR(to[j][0].z()/tozDi7[j*2],1.0,preclim);
    if ( to[j][1].z() != 0) EXPECT_NEAR(to[j][1].z()/tozDi7[j*2+1],1.0,preclim);
  }
}


TEST_F(TBDUTest, threeBD)
{
  vector<Molecule> mol;
  const int num = 7;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 5.0 ),Pt( -5.0, 0.0, 5.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < num; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 10;
  Constants const_( INTERNAL );
  const_.set_salt_concentration(0.1);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();
  
  threeBodTest.solveNmer(2, 1e-4);
  threeBodTest.solveNmer(3, 1e-4);
  threeBodTest.calcTBDEnForTor();
  
  vector<vector<double > > en = threeBodTest.getTrEn();
  vector<vector<Pt > > fr = threeBodTest.getTrFo();
  vector<vector<Pt > > to = threeBodTest.getTrTo();
  
  int j;
  for ( j = 0; j < en.size(); j++)
  {
    if ( en[j][0] != 0) EXPECT_NEAR( en[j][0]/nrgTr7[j*3], 1.0, preclim);
    if ( en[j][1] != 0) EXPECT_NEAR( en[j][1]/nrgTr7[j*3+1], 1.0, preclim);
    if ( en[j][2] != 0) EXPECT_NEAR( en[j][2]/nrgTr7[j*3+2], 1.0, preclim);
    
    if ( fr[j][0].x() != 0) EXPECT_NEAR(fr[j][0].x()/foxTr7[j*3],1.0,preclim);
    if ( fr[j][1].x() != 0) EXPECT_NEAR(fr[j][1].x()/foxTr7[j*3+1],1.0,preclim);
    if ( fr[j][2].x() != 0) EXPECT_NEAR(fr[j][2].x()/foxTr7[j*3+2],1.0,preclim);
    
    if ( fr[j][0].y() != 0) EXPECT_NEAR(fr[j][0].y()/foyTr7[j*3],1.0,preclim);
    if ( fr[j][1].y() != 0) EXPECT_NEAR(fr[j][1].y()/foyTr7[j*3+1],1.0,preclim);
    if ( fr[j][2].y() != 0) EXPECT_NEAR(fr[j][2].y()/foyTr7[j*3+2],1.0,preclim);
    
    if ( fr[j][0].z() != 0) EXPECT_NEAR(fr[j][0].z()/fozTr7[j*3],1.0,preclim);
    if ( fr[j][1].z() != 0) EXPECT_NEAR(fr[j][1].z()/fozTr7[j*3+1],1.0,preclim);
    if ( fr[j][2].z() != 0) EXPECT_NEAR(fr[j][2].z()/fozTr7[j*3+2],1.0,preclim);
    
    if ( to[j][0].x() != 0) EXPECT_NEAR(to[j][0].x()/toxTr7[j*3],1.0,preclim);
    if ( to[j][1].x() != 0) EXPECT_NEAR(to[j][1].x()/toxTr7[j*3+1],1.0,preclim);
    if ( to[j][2].x() != 0) EXPECT_NEAR(to[j][2].x()/toxTr7[j*3+2],1.0,preclim);
    
    if ( to[j][0].y() != 0) EXPECT_NEAR(to[j][0].y()/toyTr7[j*3],1.0,preclim);
    if ( to[j][1].y() != 0) EXPECT_NEAR(to[j][1].y()/toyTr7[j*3+1],1.0,preclim);
    if ( to[j][2].y() != 0) EXPECT_NEAR(to[j][2].y()/toyTr7[j*3+2],1.0,preclim);
    
    if ( to[j][0].z() != 0) EXPECT_NEAR(to[j][0].z()/tozTr7[j*3],1.0,preclim);
    if ( to[j][1].z() != 0) EXPECT_NEAR(to[j][1].z()/tozTr7[j*3+1],1.0,preclim);
    if ( to[j][2].z() != 0) EXPECT_NEAR(to[j][2].z()/tozTr7[j*3+2],1.0,preclim);
  }
}


TEST_F(TBDUTest, threeBDfor3)
{
  vector<Molecule> mol;
  const int num = 3;
  Pt pos[3] = { Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 )};
  for (int molInd = 0; molInd < num; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 10;
  Constants const_( INTERNAL );
  const_.set_salt_concentration(0.1);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  threeBodTest.solveNmer(2); threeBodTest.solveNmer(3);
  threeBodTest.calcTBDEnForTor();
  
  shared_ptr<ASolver> aSolvall = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                       make_shared<Constants>
                                                       (const_), vals);
  aSolvall->solve_A(1e-4); aSolvall->solve_gradA(1e-4);
  PhysCalc calcphys( aSolvall);
  calcphys.calc_all();
  
  int j;
  for ( j = 0; j < num; j++)
  {
    if ( calcphys.get_omegai(j) != 0)
      EXPECT_NEAR( threeBodTest.get_energyi_approx(j)/calcphys.get_omegai(j),
                  1.0, preclim);
    if ( calcphys.get_forcei(j)[0] != 0)
      EXPECT_NEAR( threeBodTest.get_forcei_approx(j).x() /
                  calcphys.get_forcei(j)[0], 1.0, preclim);
    if ( calcphys.get_forcei(j)[1] != 0)
      EXPECT_NEAR( threeBodTest.get_forcei_approx(j).y() /
                  calcphys.get_forcei(j)[1], 1.0, preclim);
    if ( calcphys.get_forcei(j)[2] != 0)
      EXPECT_NEAR( threeBodTest.get_forcei_approx(j).z() /
                  calcphys.get_forcei(j)[2], 1.0, preclim);
    
    if ( calcphys.get_taui(j)[0] != 0)
      EXPECT_NEAR( threeBodTest.get_torquei_approx(j).x() /
                  calcphys.get_taui(j)[0], 1.0, preclim);
    if ( calcphys.get_taui(j)[1] != 0)
      EXPECT_NEAR( threeBodTest.get_torquei_approx(j).y() /
                  calcphys.get_taui(j)[1], 1.0, preclim);
    if ( calcphys.get_taui(j)[2] != 0)
      EXPECT_NEAR( threeBodTest.get_torquei_approx(j).z() /
                  calcphys.get_taui(j)[2], 1.0, preclim);
  }
}

TEST_F(TBDUTest, threeBDfor7)
{
  vector<Molecule> mol;
  const int num = 7;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < num; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 10;
  Constants const_( INTERNAL );
  const_.set_salt_concentration(0.1);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();
  
  threeBodTest.solveNmer(2, 1e-5);
  threeBodTest.solveNmer(3, 1e-5);
  threeBodTest.calcTBDEnForTor();
  
  int j;
  for ( j = 0; j < num; j++)
  {
    if ( threeBodTest.get_energyi_approx(j) != 0)
      EXPECT_NEAR( threeBodTest.get_energyi_approx(j)/nrg7[j], 1.0, preclim);
    if ( threeBodTest.get_forcei_approx(j).x() != 0)
      EXPECT_NEAR( threeBodTest.get_forcei_approx(j).x()/fox7[j], 1.0, preclim);
    if ( threeBodTest.get_forcei_approx(j).y() != 0)
      EXPECT_NEAR( threeBodTest.get_forcei_approx(j).y()/foy7[j], 1.0, preclim);
    if ( threeBodTest.get_forcei_approx(j).z() != 0)
      EXPECT_NEAR( threeBodTest.get_forcei_approx(j).z()/foz7[j], 1.0, preclim);
    
    if ( threeBodTest.get_torquei_approx(j).x() != 0)
      EXPECT_NEAR( threeBodTest.get_torquei_approx(j).x()/tox7[j],1.0,preclim);
    if ( threeBodTest.get_torquei_approx(j).y() != 0)
      EXPECT_NEAR( threeBodTest.get_torquei_approx(j).y()/toy7[j],1.0,preclim);
    if ( threeBodTest.get_torquei_approx(j).z() != 0)
      EXPECT_NEAR( threeBodTest.get_torquei_approx(j).z()/toz7[j],1.0,preclim);
  }
  cout << endl;
}

#endif /* ThreeBodyUnitTest_h */
