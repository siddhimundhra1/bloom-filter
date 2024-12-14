
#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <fstream>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <map>
using namespace std;
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;



class Bloomfilter{
private:
    int N;
    int m;
    int k;
    unsigned int p=2147483647;
    mt19937 mt;
    uniform_int_distribution<int> dist_a;
    uniform_int_distribution<int> dist_b;
    vector<int> a, b, s;
    vector<bool> filter;


public:
    Bloomfilter(int N, int m, int k)
        : N(N), m(m), k(k), dist_a(1, p - 1), dist_b(0, p - 1){
        a.resize(k);
        b.resize(k);
        s.resize(k);
        filter.resize(m);
        mt.seed(static_cast<unsigned int>(time(nullptr)));

        for (int i = 0; i < k; ++i) {
            a[i] = dist_a(mt);
            b[i] = dist_b(mt);
            s[i]=time(nullptr)+i;
        }
    }
    int hashfunction1(int x, int i) {
        return ((a[i]*x+b[i])%p)%m;
    }

    int hashfunction2(int x, int i)
    {

        mt19937 rng(s[i]+x);
        uniform_int_distribution<int> dist(0, 2147483647);
        return dist(rng) % m;
    }
    void add1(int x)
    {
        for (int i=0; i<k; ++i) {
            int hash_value = hashfunction1(x,i);
            filter[hash_value]=true;
        }

    }
    bool contains1(int x)
    {
        for (int i = 0; i < k; ++i) {
            int hash_value=hashfunction1(x, i);
            if (!filter[hash_value]) {
                return false;
            }
        }
        return true;
    }
    void add2(int x)
    {
        for (int i=0; i<k; ++i) {
            int hash_value = hashfunction2(x,i);
            filter[hash_value]=true;
        }

    }
    bool contains2(int x)
    {
        for (int i = 0; i < k; ++i) {
            int hash_value=hashfunction2(x, i);
            if (!filter[hash_value]) {
                return false;
            }
        }
        return true;
    }

    


};

double testFalsePositiveRate1(Bloomfilter& bf, const unordered_set<int>& inserted, int n_queries, int N) {
    mt19937 mt(static_cast<unsigned int>(time(nullptr)));
    uniform_int_distribution<int> universe(0, N - 1);

    int false_positives = 0;
    for (int i=0; i<n_queries;++i) {
        int query = universe(mt);
        if (inserted.find(query) == inserted.end() && bf.contains1(query)) {
            false_positives++;
        }
    }
    return static_cast<double>(false_positives)/n_queries;

}


double testFalsePositiveRate2(Bloomfilter& bf, const unordered_set<int>& inserted, int n_queries, int N) {
    mt19937 mt(static_cast<unsigned int>(time(nullptr)));
    uniform_int_distribution<int> universe(0,N);

    int false_positives = 0;
    for (int i=0;i<n_queries;++i) {
        int query = universe(mt);
        if (inserted.find(query) == inserted.end() && bf.contains2(query)) {
            false_positives++;
        }
    }
    return static_cast<double>(false_positives) / n_queries;

}


int main() {

    int N = 1073741824;      
    int n=1000000;
    int k=1;
    int m = 100000;
    Bloomfilter bf(N, m, k);

    

    // Random data
    mt19937 mt;
    uniform_int_distribution<int> universe(0, N - 1);

    vector<int> collisions_per_bucket1(m,0);
    vector<int> collisions_per_bucket2(m,0);
    for (int i=0; i<n; ++i)//comment this out when using seq
    {
        int val=universe(mt);
        int m_i_1=bf.hashfunction1(val, 0); //gives bucket #
        int m_i_2=bf.hashfunction2(val, 0);
        collisions_per_bucket1[m_i_1]+=1;
        collisions_per_bucket2[m_i_2]+=1;

    }
    /*for (int i=0; i<2*n; i=i+2) //sequential data with correlation
    {
        
        int m_i_1=bf.hashfunction1(i, 0); //gives bucket #
        int m_i_2=bf.hashfunction2(i, 0);
        collisions_per_bucket1[m_i_1]+=1;
        collisions_per_bucket2[m_i_2]+=1;

    }*/

    unordered_map<int,int> collisionfreq1; // # of collisions vs freq
    unordered_map<int,int> collisionfreq2;

    for (int i: collisions_per_bucket1)
    {
        if (collisionfreq1.find(i)==collisionfreq1.end())
            collisionfreq1[i]=1;
        else
            collisionfreq1[i]+=1;
    }

    for (int i: collisions_per_bucket2)
    {
        if (collisionfreq2.find(i)==collisionfreq2.end())
            collisionfreq2[i]=1;
        else
            collisionfreq2[i]+=1;
    }
    vector<pair<int, int>> sorted_collisionfreq1(collisionfreq1.begin(), collisionfreq1.end());
    vector<pair<int, int>> sorted_collisionfreq2(collisionfreq2.begin(), collisionfreq2.end());

    sort(sorted_collisionfreq1.begin(), sorted_collisionfreq1.end());
    sort(sorted_collisionfreq2.begin(), sorted_collisionfreq2.end());

    vector<int> xcollisions1, xcollisions2, yfrequency1, yfrequency2;
    for (auto p: sorted_collisionfreq1)
    {
        xcollisions1.push_back(p.first);
        yfrequency1.push_back(p.second);
    }

    for (auto p: sorted_collisionfreq2)
    {
        xcollisions2.push_back(p.first);
        yfrequency2.push_back(p.second);
    }

    plt::figure();
    plt::plot(xcollisions1, yfrequency1, {{"label", "Hash Function 1"}, {"color", "blue"}});
    plt::plot(xcollisions2, yfrequency2, {{"label", "Hash Function 2"}, {"color", "red"}});
    plt::title("Number of Collisions vs Number of Buckets for random data");
    plt::xlabel("Number of Collisions");
    plt::ylabel("Number of Buckets");
    plt::legend();
    plt::show();

    return 0;
    
    
    /*int N = 1073741824;  
    int c = 15;     
    int n=100000;
    int m = c*n;  
    int n_queries = 100000; 


    vector<int> k_values = {4,5,6,7,8,9,10,11,12,13,14};
    sort(k_values.begin(),k_values.end());
    vector<double> false_positive_rates1;
    vector<double> false_positive_rates2;

    mt19937 mt(static_cast<unsigned int>(time(nullptr)));
    uniform_int_distribution<int> universe(0, N);
    for (int k : k_values) {
        vector<double> trial_results;

        for (int trial = 0; trial < 10; ++trial) {
            Bloomfilter bf(N, m, k);

            unordered_set<int> inserted_elements;
            for (int i = 0; i < n; ++i) {
                int query = universe(mt);
                bf.add1(query);
                inserted_elements.insert(query);
            }

            
            double fpr = testFalsePositiveRate1(bf, inserted_elements, n_queries, N);
            trial_results.push_back(fpr);
        }

        
        sort(trial_results.begin(), trial_results.end());
        double median_fpr = trial_results[trial_results.size() / 2];
        false_positive_rates1.push_back(median_fpr);

        cout << "k = " << k << ", Median False Positive Rate HF#1 = " << median_fpr << " Theoretical: " << pow(1 - exp(-1.0*k*n/m),k) <<endl;
    }



    for (int k : k_values) {
        vector<double> trial_results;

        for (int trial = 0; trial < 10; ++trial) {
            Bloomfilter bf(N, m, k);

            unordered_set<int> inserted_elements;
            for (int i = 0; i < n; ++i) {
                int query = universe(mt);
                bf.add2(query);
                inserted_elements.insert(query);
            }

            double fpr = testFalsePositiveRate2(bf, inserted_elements, n_queries, N);
            trial_results.push_back(fpr);
        }

        sort(trial_results.begin(), trial_results.end());
        double median_fpr = trial_results[trial_results.size() / 2];
        false_positive_rates2.push_back(median_fpr);

        cout << "k = " << k << ", Median False Positive Rate HF#2 = " << median_fpr << " Theoretical: " << pow(1 - exp(-1.0 * k * n / m),k) <<endl;
    }



vector<double> theoretical_fpr;
for (int k : k_values) {
    double fpr = pow(1-exp(-1.0*k*n/m),k);
    theoretical_fpr.push_back(fpr);
}

    plt::figure();
    plt::named_plot("Experimental Rate HF#1", k_values, false_positive_rates1, "o-"); // Experimental Hash function #1
    plt::named_plot("Experimental Rate HF#2", k_values, false_positive_rates2, "o-"); // Experimental Hash function #2
    plt::named_plot("Theoretical Rate", k_values, theoretical_fpr, "x-");      // Theoretical FPR

    
    plt::xlabel("Number of Hash Functions (k)");
    plt::ylabel("False Positive Rate");
    plt::title("False Positive Rate vs Number of Hash Functions for c=15");

    plt::legend(); 
    plt::show();

    return 0;*/
}
