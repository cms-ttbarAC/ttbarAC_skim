/*
Created:        14 May       2016
Last Updated:   20 August    2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Common tools needed

*/
#include "Analysis/CyMiniAna/interface/tools.h"

namespace cma{


void check_file(const std::string & filename) {
    /* Check that file exists */
    std::ifstream f = open_file(filename);

    return;
}


void check_file(const std::ifstream& file, const std::string& fname){
    /* Check that file exists */

    if (!file) {
        cma::ERROR("TOOLS : File does not exist:       "+fname);
        cma::ERROR("TOOLS : Exiting. ");
        assert(file);
    }

    return;
}

std::ifstream open_file(const std::string &filename) {
    /* Open file */
    std::ifstream ifile(filename.c_str());
    check_file(ifile,filename);

    return ifile;
}


void read_file( const std::string &file_name, std::vector<std::string> &values, const std::string &comment ) {
    /* Read in a generic file and put it into a vector of strings */
    std::ifstream tmp_name = open_file(file_name);

    // open the file and put the data into a vector
    std::string line("");
    if (tmp_name.is_open()){

        while (std::getline(tmp_name, line)) {
            std::string newstring(line);

            // allow for comments
            std::size_t lineComment = line.find(comment);
            if (lineComment != std::string::npos) newstring = line.substr(0, lineComment);

            // remove all white spaces at the end of the string
            std::size_t space_pos = newstring.rfind(" ");
            while ( space_pos != std::string::npos && space_pos == newstring.size()-1 ) {
                newstring = newstring.substr(0, newstring.rfind(" "));
                space_pos = newstring.rfind(" ");
            }

            // ignore empty lines
            if(newstring.length()==0) continue;

            values.push_back(newstring); // put values into vector
        }

        tmp_name.close();
    }

    return;
}



void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }

    return;
}


void getListOfKeys( TFile* file, std::vector<std::string> &fileKeys ){
    /* Find the list of TTrees in a file */
    fileKeys.clear();

    TList *list = file->GetListOfKeys();
    TIter iter(list->MakeIterator());
    while(TObject* obj = iter()){
        TKey* key = (TKey*)obj;
        std::string keyname( key->GetName() );
        fileKeys.push_back(keyname);
    }

    return;
}


void getSampleWeights( std::string metadata_file,
                       std::map<std::string,Sample>& samples){
    /* Calculate XSection, KFactor, NEvents, and sum of weights (AMI) */
    cma::INFO("TOOLS : Get sample weights (including sum of weights)");

    // get the absolute path (in case of batch job)
    char* cma_path = getenv("CYMINIANADIR");
    std::string cma_absPath("");
    if (cma_path==NULL){
        cma::WARNING("TOOLS : environment variable 'CYMINIANADIR' is not set." );
        cma::WARNING("TOOLS : Relative paths will be used " );
        cma_absPath = "./";
    }
    else cma_absPath = cma_path;


    std::ifstream in( (cma_absPath+"/"+metadata_file).c_str());
    if (!in) cma::WARNING("TOOLS : File does not exist: "+cma_absPath+"/"+metadata_file);

    std::string line;
    samples.clear();
    while( std::getline(in,line) ) {

        if (!line.empty() && line[0]!='#') {
            std::string dsid("");
            unsigned int NEvents;
            float xSect,kFact,sumWeights;

            std::istringstream istr(line);
            istr >> dsid >> xSect >> sumWeights >> kFact >> NEvents;

            Sample s;
            s.primaryDataset = dsid;
            s.XSection       = xSect;
            s.KFactor        = kFact;
            s.NEvents        = NEvents;
            s.sumOfWeights   = sumWeights;

            samples[dsid] = s;
        }
    }
    in.close();

    return;
}


bool str2bool( const std::string value ){
    /* Turn string into boolean */
    bool valueBoolean(false);

    if (value.compare("True")==0 || value.compare("true")==0 || value.compare("1")==0){
        valueBoolean = true;
    }
    else{
        valueBoolean = false;
    }

    return valueBoolean;
}


std::string vectorToStr( const std::vector<std::string> &vec ){
    std::string str_list;
    for( const std::string &str : vec)
        str_list += str + std::string(",");
    str_list.pop_back(); // remove last comma
    return str_list;
}


unsigned int setRandomNumberSeeds(const Lepton& lepton, const Lepton& antiLepton, 
                                  const Jet& jet1, const Jet& jet2) {
    /* 
       Asymmetric treatment of both jets, and also both leptons, 
       to ensure different seed for each combination in dileptonTtbarReco
    */
    unsigned int seed = static_cast<int>( 1.e6 * (jet1.p4.Pt()/jet2.p4.Pt()) * 
                                          std::sin((lepton.p4.Pt() + 2.*antiLepton.p4.Pt()) * 1.e6) );
    gRandom->SetSeed(seed);

    return seed;
}


bool deltaRMatch( const TLorentzVector &particle1, const TLorentzVector &particle2, const double deltaR ){
    /* Do the deltaR calculation (in one place) */
    return (particle1.DeltaR(particle2)<deltaR);
}


std::string m_debugLevel = "SetMe";
void setVerboseLevel( const std::string& verboseLevel ) {
    m_debugLevel = verboseLevel;
    return;
}

void DEBUG(const std::string& message){
    /* Debug level (verbosity of output) */
    verbose("DEBUG",message);
    return;
}
void INFO(const std::string& message){
    /* Info level (verbosity of output) */
    verbose("INFO",message);
    return;
}
void WARNING(const std::string& message){
    /* Warning level (verbosity of output) */
    verbose("WARNING",message);
    return;
}
void ERROR(const std::string& message){
    /* Error level (verbosity of output) */
    verbose("ERROR",message);
    return;
}

void verbose(const std::string level, const std::string& message){
    /* 
       Printing output to console (debug,warning,error messages)
         if the level is "DEBUG", then all messages should be printed (DEBUG/INFO/WARNING/ERROR)
         if the level is "INFO", then only INFO/WARNING/ERROR messages should be printed
         if the level is "WARNING", then only WARNING/ERROR messages should be printed
         if the level is "ERROR", then only ERROR messages should be printed
    */
    std::map<std::string,unsigned int> debugMap = {
            {"DEBUG",   0},
            {"INFO",    1},
            {"WARNING", 2},
            {"ERROR",   3} };

    if ( debugMap.at( level ) >= debugMap.at( m_debugLevel ))
        std::cout << " " << level << " :: " << message << std::endl;

    return;
}

std::map<std::string,unsigned int> verboseMap() {
    /* mapping of verbose level to integer */
    std::map<std::string,unsigned int> verbose_map = {
            {"DEBUG",   0},
            {"INFO",    1},
            {"WARNING", 2},
            {"ERROR",   3} };
    
    return verbose_map;
}

void HELP(const std::string& runExecutable){
    /* HELP message (pass 'runExecutable' in case you are running from some 
       script like 'skim', 'run', or a custom macro)
    */
    std::cout << "\n   ** CyMiniAna ** " << std::endl;
    std::cout << "   --------------- " << std::endl;
    std::cout << "   Framework to perform event selection, write-out" << std::endl;
    std::cout << "   a few histograms or efficiencies, and make plots.\n" << std::endl;

    std::cout << "   To run:" << std::endl;
    std::cout << "      ./" << runExecutable << " share/cmaConfig.txt \n" << std::endl;
    std::cout << "    where 'share/cmaConfig.txt' is the configuration file \n" << std::endl;

    return;
}

} // end namespace

// the end

