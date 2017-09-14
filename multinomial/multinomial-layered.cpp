// Code shared with IsoSpec
//
//

#include <cmath>
#include <vector>
#include <cstdlib>
#include <limits>
#include <cstring>
#include <unordered_set>
#include <algorithm>




typedef int* Conf;


static inline double logFactorial(int n) { return lgamma(n+1); }

inline double logProb(const int* conf, const double* logProbs, int dim)
{
    int     N = 0;
    double  res = 0.0;

    for(int i=0; i < dim; i++)
    {
        N   += conf[i];
        res -= logFactorial(conf[i]);
        res += conf[i] * logProbs[i];
    }
    return res + logFactorial(N);
}



class Marginal
{
private:
    bool disowned;
protected:
    const unsigned int isotopeNo;
    const unsigned int atomCnt;
    const double* atom_lProbs;
    const Conf mode_conf;
    const double mode_lprob;

public:
    Marginal(
        const double* _probs,
        int _isotopeNo,                  // No of isotope configurations.
        int _atomCnt
    );
    Marginal(Marginal& other) = delete;
    Marginal& operator= (const Marginal& other) = delete;
    Marginal(Marginal&& other);
    virtual ~Marginal();

    inline int get_isotopeNo() const { return isotopeNo; };
    inline double getModeLProb() const { return mode_lprob; };
};



Conf initialConfigure(const int atomCnt, const int isotopeNo, const double* probs, const double* lprobs)
{
    Conf res = new int[isotopeNo];

    for(int i = 0; i < isotopeNo; ++i )
    {
        res[i] = int( atomCnt * probs[i] ) + 1;
    }

    int s = 0;

    for(int i = 0; i < isotopeNo; ++i) s += res[i];

    int diff = atomCnt - s;

    // Too little: enlarging fist index.
    if( diff > 0 ){
        res[0] += diff;
    }
    // Too much: trying to redistribute the difference: hopefully the first element is the largest.
    if( diff < 0 ){
        diff = abs(diff);
        int i = 0, coordDiff = 0;

        while( diff > 0){
            coordDiff = res[i] - diff;

            if( coordDiff >= 0 ){
                res[i] -= diff;
                diff = 0;
            } else {
                res[i] = 0;
                i++;
                diff = abs(coordDiff);
            }
        }
    }

    // What we computed so far will be very close to the mode: hillclimb the rest of the way

    bool modified = true;
    double LP = logProb(res, lprobs, isotopeNo);
    double NLP;

    while(modified)
    {
    	modified = false;
    	for(int ii = 0; ii<isotopeNo; ii++)
	    for(int jj = 0; jj<isotopeNo; jj++)
	        if(ii != jj and res[ii] > 0)
		{
		    res[ii]--;
		    res[jj]++;
		    NLP = logProb(res, lprobs, isotopeNo);
		    if(NLP>LP or (NLP==LP and ii>jj))
		    {
		    	modified = true;
			LP = NLP;
		    }
		    else
		    {
		    	res[ii]++;
			res[jj]--;
		    }
		}


    }
    return res;
}



template <typename T> inline static T* array_copy(const T* A, int size)
{
    T* ret = new T[size];
    memcpy(ret, A, size*sizeof(T));
    return ret;
}

double* getMLogProbs(const double* probs, int isotopeNo)
{
    double* ret = new double[isotopeNo];
    for(int ii=0; ii<isotopeNo; ii++)
        ret[ii] = log(probs[ii]);
    return ret;
}

template <typename T> inline void copyConf(
    const T* source, T* destination,
    int dim
){
    memcpy(destination, source, dim*sizeof(T));
}

template <typename T> class Allocator{
private:
    T*      currentTab;
    int currentId;
    const int       dim, tabSize;
    std::vector<T*>  prevTabs;
public:
    Allocator(const int dim, const int tabSize = 10000);
    ~Allocator();

    void shiftTables();

    inline T* newConf()
    {
        currentId++;

        if (currentId >= tabSize)
            shiftTables();

        return &(currentTab[ currentId * dim ]);
    }

    inline T* makeCopy(const T* conf)
    {
        T* currentPlace = newConf();
        copyConf<T>( conf, currentPlace, dim );

        return currentPlace;
    }

    inline T* makeExternalCopy(const T* conf)
    {
        T* res = new T[dim];
        copyConf( conf, res, dim );

        return res;
    }
};

template <typename T>
Allocator<T>::Allocator(const int dim, const int tabSize): currentId(-1), dim(dim), tabSize(tabSize)
{
    currentTab = new T[dim * tabSize];
}

template <typename T>
Allocator<T>::~Allocator()
{
    for(unsigned int i = 0; i < prevTabs.size(); ++i)
    {
        delete [] prevTabs[i];
    }

    delete [] currentTab;
}

template <typename T>
void Allocator<T>::shiftTables()
{
    prevTabs.push_back(currentTab);
    currentTab      = new T[dim * tabSize];
    currentId       = 0;
}
template class Allocator<int>;


class KeyHasher
{
private:
    int dim;
public:
    KeyHasher(int dim);

    inline std::size_t operator()(const int* conf) const
    {
        // Following Boost...
        std::size_t seed = 0;
        for(int i = 0; i < dim; ++i )
            seed ^= conf[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    };
};


class ConfEqual
{
private:
    int size;
public:
    ConfEqual(int dim);

    inline bool operator()(const int* conf1, const int* conf2) const
    {
        return not memcmp(conf1, conf2, size);
    }
};

class ConfOrderMarginalDescending
{
//configurations comparator
    const double*  logProbs;
    int dim;
public:
    ConfOrderMarginalDescending(const double* logProbs, int dim);

    inline bool operator()(const Conf conf1, const Conf conf2)
    {// Return true if conf1 is less probable than conf2.
        return logProb(conf1,logProbs,dim) > logProb(conf2,logProbs,dim);
    };
};


KeyHasher::KeyHasher(int dim)
: dim(dim)
{}

ConfEqual::ConfEqual(int dim)
: size( dim*sizeof(int) )
{}

ConfOrderMarginalDescending::ConfOrderMarginalDescending(const double* logProbs, int dim)
: logProbs(logProbs), dim(dim)
{}



Marginal::Marginal(
    const double* _probs,
    int _isotopeNo,
    int _atomCnt
) :
disowned(false),
isotopeNo(_isotopeNo),
atomCnt(_atomCnt),
atom_lProbs(getMLogProbs(_probs, isotopeNo)),
mode_conf(initialConfigure(atomCnt, isotopeNo, _probs, atom_lProbs)),
mode_lprob(logProb(mode_conf, atom_lProbs, isotopeNo))
{}

Marginal::Marginal(Marginal&& other) :
disowned(other.disowned),
isotopeNo(other.isotopeNo),
atomCnt(other.atomCnt),
atom_lProbs(other.atom_lProbs),
mode_conf(other.mode_conf),
mode_lprob(logProb(mode_conf, atom_lProbs, isotopeNo))
{
    other.disowned = true;
}

Marginal::~Marginal()
{
    if(not disowned)
    {
        delete[] atom_lProbs;
        delete[] mode_conf;
    }
}



class LayeredMarginal : public Marginal
{
    double current_threshold, new_threshold;
    std::vector<Conf> configurations;
    std::vector<Conf> fringe;
    Allocator<int> allocator;
    unsigned int sorted_up_to_idx;
    const ConfEqual equalizer;
    const KeyHasher keyHasher;
    const ConfOrderMarginalDescending orderMarginal;
    const int hashSize;
    std::unordered_set<Conf,KeyHasher,ConfEqual> visited;
    double opc;
    Conf currentConf;
    std::vector<Conf> new_fringe;



public:
    LayeredMarginal(Marginal&& m, int tabSize = 1000, int hashSize = 1000);
    bool extend(double new_threshold);
    bool next();
    inline double get_lProb() const { return opc; };
    inline const Conf& get_conf() const { return currentConf; };

};





LayeredMarginal::LayeredMarginal(Marginal&& m, int tabSize, int _hashSize)
: Marginal(std::move(m)), current_threshold(1.0), allocator(tabSize), sorted_up_to_idx(0),
equalizer(isotopeNo), keyHasher(isotopeNo), orderMarginal(atom_lProbs, isotopeNo), hashSize(_hashSize),
visited(hashSize,keyHasher,equalizer)
{
    fringe.push_back(mode_conf);
}

bool LayeredMarginal::extend(double _new_threshold)
{
    if(fringe.empty())
        return false;

    new_threshold = _new_threshold;

    // TODO: Make sorting optional (controlled by argument?)
    new_fringe.clear();

    visited.clear();

    for(unsigned int ii = 0; ii<fringe.size(); ii++)
        visited.insert(fringe[ii]);

    current_threshold = new_threshold;
    fringe.swap(new_fringe);

    return true;
}

bool LayeredMarginal::next()
{
    double lpc;

    while(true)
    {
        if(fringe.empty())
            return false;

        currentConf = fringe.back();
        fringe.pop_back();

        opc = logProb(currentConf, atom_lProbs, isotopeNo);
        if(opc >= new_threshold)
            configurations.push_back(currentConf);
        else
            break;
    }


    for(unsigned int ii = 0; ii < isotopeNo; ii++ )
        for(unsigned int jj = 0; jj < isotopeNo; jj++ )
            if( ii != jj and currentConf[jj] > 0 )
            {
                currentConf[ii]++;
                currentConf[jj]--;

                lpc = logProb(currentConf, atom_lProbs, isotopeNo);

                if (visited.count(currentConf) == 0 and lpc < current_threshold and
                    (opc > lpc or (opc == lpc and ii > jj)))
                {
                    Conf nc = allocator.makeCopy(currentConf);
                    visited.insert(nc);
                    if(lpc >= new_threshold)
                        fringe.push_back(nc);
                    else
                        new_fringe.push_back(nc);
                }

                currentConf[ii]--;
                currentConf[jj]++;

            }

    return true;
}



