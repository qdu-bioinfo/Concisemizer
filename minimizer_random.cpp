// minimizer_random sampling
// compile：g++ -std=c++11 minimizer_random.cpp -o minimizer_random
// operation：./minimizer_random -n XXX.fna -k xxx -w xxx
// "-n XXX.fna" reference genome file
// "-k xxx" k-mer length
// "-w xxx" window size
// Output five lines of results. The first line is k and w; Second row density; The third row density factor; The fourth line: E-hits; The fifth line: DR-Score.

#include <iostream>
#include <limits>  // Involving numeric_limits, the extreme values of various variable types.
#include <utility> // pair
#include <fstream>
#include <string>
#include <deque> // Insert at the end of the queue, delete from the front of the queue, traverse
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
using namespace std;

static inline unsigned long long int hash64(unsigned long long int key, unsigned long long int mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}
static inline unsigned long long int hash64_2(unsigned long long int key, unsigned long long int mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)); // Since the lexicographical order also needs to be made into a hash table, the &mask was removed to facilitate comparison with the lexicographical order.
    return key;
}

int main(int argc, char *argv[])
{
    int k_mer_length = 0; // Kmer length.
    int k_mer_win = 0;    // Window size.

    // argc represents the number of parameters, including the program name.
    // argv is an array of strings, containing each parameter.
    if (argc != 7)
    {
        std::cout << "Error: No arguments provided." << std::endl
                  << "Please run the program in the following format:" << std::endl
                  << "./minimizer_random -n XXX.fna -k xxx -w xxx" << std::endl
                  << std::endl;
        return 0;
    }
    // Obtain option
    std::string option_n = argv[1];
    std::ifstream file(argv[2]);
    // Evaluate the option
    if (option_n == "-n")
    {
        if (!file)
        {
            // If the file fails to open, print the error message and return.
            std::cerr << "Failed to open file: '" << argv[2] << "', please check whether the file is damaged." << std::endl;
            return 0;
        }
        // Obtain option
        std::string option_k = argv[3];
        // Evaluate the option
        if (option_k == "-k")
        {
            try
            {
                k_mer_length = stoi(argv[4]); // Recommended usage: Automatically check for illegal input.
            }
            catch (exception &e)
            {
                cerr << "Invalid number: " << e.what() << endl;
                return 0;
            }
            // Obtain option
            std::string option_w = argv[5];
            // Evaluate the option
            if (option_w == "-w")
            {
                try
                {
                    k_mer_win = stoi(argv[6]); // Recommended usage: Automatically check for illegal input.
                }
                catch (exception &e)
                {
                    cerr << "Invalid number: " << e.what() << endl;
                    return 0;
                }
            }
            else
            {
                std::cout << "Error: No arguments provided." << std::endl
                          << "Please run the program in the following format:" << std::endl
                          << "./minimizer_random -n XXX.fna -k xxx -w xxx" << std::endl
                          << std::endl;
                return 0;
            }
        }
        else
        {
            std::cout << "Error: No arguments provided." << std::endl
                      << "Please run the program in the following format:" << std::endl
                      << "./minimizer_random -n XXX.fna -k xxx -w xxx" << std::endl
                      << std::endl;
            return 0;
        }
    }
    else
    {
        std::cout << "Error: No arguments provided." << std::endl
                  << "Please run the program in the following format:" << std::endl
                  << "./minimizer_random -n XXX.fna -k xxx -w xxx" << std::endl
                  << std::endl;
        return 0;
    }
    // To prevent kmer and its transposition from being the same, set the length of kmer to an odd number.
    if (k_mer_length < 5 || k_mer_length > 31)
    {
        std::cout << "Please re-enter the value of k. Please note that the value of k should be within the range of [5, 31]." << std::endl
                  << std::endl;
        return 0;
    }
    if (k_mer_win < 2 || k_mer_win > 100)
    {
        std::cout << "Please re-enter the value of w. Note that the magnitude of the w value should be an integer within the range of [2,100]." << std::endl
                  << std::endl;
        return 0;
    }

    unsigned long long int ATUCG[256] = {0};
    for (int i = 0; i < 256; i++)
    {
        ATUCG[i] = 4;
    }
    ATUCG['A'] = 0;
    ATUCG['a'] = 0;
    ATUCG['C'] = 1;
    ATUCG['c'] = 1;
    ATUCG['G'] = 2;
    ATUCG['g'] = 2;
    ATUCG['T'] = 3;
    ATUCG['t'] = 3;
    ATUCG['U'] = 3;
    ATUCG['u'] = 3;

    unsigned long long int ATUCG0[256] = {0};
    for (int i = 0; i < 256; i++)
    {
        ATUCG0[i] = 4;
    }
    // The pairing relationship is A-TU and C-G, which are reverse complementary sequences. When complementing, simply invert them directly.
    ATUCG0['A'] = (((unsigned long long int)3) << 62);
    ATUCG0['a'] = (((unsigned long long int)3) << 62);
    ATUCG0['C'] = (((unsigned long long int)2) << 62);
    ATUCG0['c'] = (((unsigned long long int)2) << 62);
    ATUCG0['G'] = (((unsigned long long int)1) << 62);
    ATUCG0['g'] = (((unsigned long long int)1) << 62);
    ATUCG0['T'] = (((unsigned long long int)0) << 62);
    ATUCG0['t'] = (((unsigned long long int)0) << 62);
    ATUCG0['U'] = (((unsigned long long int)0) << 62);
    ATUCG0['u'] = (((unsigned long long int)0) << 62);

    pair<unsigned long long int, unsigned long long int> minimizer_mini; // The value of mini.
    minimizer_mini.first = numeric_limits<unsigned long long int>::max();
    minimizer_mini.second = numeric_limits<unsigned long long int>::max();
    unsigned long long int minimizer_minihash; // The hash value of mini
    minimizer_minihash = numeric_limits<unsigned long long int>::max();
    int dq_min_bj = 1; // The position of the minimum element in the deque, random index access flag

    unsigned long long int mask = (1ULL << 2 * k_mer_length) - 1;
    if (k_mer_length == 32)
    {
        mask = 0ULL - 1;
    }

    std::cout << k_mer_length << " " << k_mer_win << std::endl;
    int weiyi = 64 - k_mer_length * 2;
    int k_mer_bj = k_mer_length; // Mark whether kmer is available.
    string genome_line;
    unsigned long long int x = 0;  // Kmer value
    unsigned long long int x0 = 0; // Complementary pairing of kmer
    unsigned long long int x3 = 0; // Store the maximum and minimum values in x and x0.
    unsigned long long int line_n = 0;
    unsigned long long int suoyin = 0;     // The range of the index depends on the size of the largest reference genome.
                                           // The starting position of the index may depend on the maximum length that the sequencing instrument can read.
                                           // During the alignment process, the size can also be limited to prevent subtraction from resulting in a negative number.
    unsigned long long int s_kmer_sum = 0; // Counting the number of k-mers in the reference sequence.

    vector<pair<unsigned long long int, unsigned long long int>> save_k_mer;  // Save the found minimizer and its location information.
    deque<pair<unsigned long long int, unsigned long long int>> dq_minimizer; // Insert at the end of the queue, delete at the beginning of the queue, and traverse.

    // Read the file, process one line at a time.
    while (getline(file, genome_line))
    {
        if (genome_line[0] != '>') // It is not the name of a microbial sequence.Entering this statement indicates that a complete kmer has already been obtained.
        {
            line_n = genome_line.length(); // The length of the string.
            for (unsigned long long int i = 0; i < line_n; i++)
            {
                suoyin++;
                if (ATUCG[genome_line[i]] != 4)
                {
                    x0 = (x0 >> 2) | ATUCG0[genome_line[i]];
                    x = (x << 2) | ATUCG[genome_line[i]];
                    x = (x << weiyi) >> weiyi;
                    s_kmer_sum++;
                    x3 = x0 >> weiyi;
                    if (x3 == x)
                    {
                        continue;
                    }
                    if (x3 > x)
                    {
                        x3 = x;
                    }
                    // Delete the earliest kmer that joined.
                    dq_minimizer.pop_front();
                    dq_min_bj--;        // Move the marker to the next one
                    if (dq_min_bj == 0) // What is deleted is the minimum value. A new minimum value must be found and saved.
                    {
                        // insert

                        dq_minimizer.push_back(make_pair(x3, suoyin));
                        // Find a new minimum value.
                        auto it = dq_minimizer.begin();
                        minimizer_mini = *it;
                        minimizer_minihash = hash64((*it).first, mask);
                        dq_min_bj = 1;
                        int i = 1;
                        while (++it != dq_minimizer.end())
                        {
                            i++;
                            if (hash64((*it).first, mask) <= minimizer_minihash) // Same as min select new.
                            {
                                minimizer_mini = *it;
                                minimizer_minihash = hash64((*it).first, mask);
                                dq_min_bj = i; // Point to the new minimum value in dp
                            }
                        }
                        save_k_mer.push_back(make_pair(hash64(minimizer_mini.first, mask), minimizer_mini.second));
                    }
                    else // What is deleted is not the minimum value. If the newly inserted element is min, it will be saved.
                    {

                        dq_minimizer.push_back(make_pair(x3, suoyin));
                        // Compare the new minimum value with the previous minimum value
                        // If they are the same, do nothing
                        // If they are different, save
                        if (hash64(x3, mask) <= minimizer_minihash) // Same as min select new.
                        {
                            minimizer_mini.first = x3;
                            minimizer_mini.second = suoyin;
                            minimizer_minihash = hash64(x3, mask);
                            dq_min_bj = k_mer_win; // Point to the new minimum value in dp
                            save_k_mer.push_back(make_pair(hash64(minimizer_mini.first, mask), minimizer_mini.second));
                        }
                    }
                }
                else // !atcgu Base, Find the new kmer deque.
                {
                    dq_minimizer.clear();
                    char ch_atcg;
                    dq_min_bj = 1;
                    minimizer_mini.first = numeric_limits<unsigned long long int>::max();
                    minimizer_mini.second = numeric_limits<unsigned long long int>::max();
                    minimizer_minihash = numeric_limits<unsigned long long int>::max();
                    x = 0;
                    x0 = 0;
                    k_mer_bj = k_mer_length;
                    while (k_mer_bj != 0)
                    {
                        if ((i + 1) >= line_n) // I have finished reading this.
                        {
                            if (file >> ch_atcg) // Have you reached the end of the file?
                            {
                                ;
                            }
                            else
                            {
                                break;
                            }
                        }
                        else
                        {
                            i++; // After that, proceed to use
                            ch_atcg = genome_line[i];
                        }
                        if (ch_atcg == '>')
                        {
                            file.putback(ch_atcg); // Backspace one character
                            break;
                        }
                        if (ch_atcg != '\n' && ch_atcg != '\r')
                        {
                            suoyin++;
                            if (ATUCG[ch_atcg] != 4)
                            {
                                x0 = (x0 >> 2) | ATUCG0[ch_atcg];
                                x = (x << 2) | ATUCG[ch_atcg];
                                x = (x << weiyi) >> weiyi;
                                k_mer_bj--;
                                if (k_mer_bj == 0)
                                {
                                    s_kmer_sum++;
                                    x3 = x0 >> weiyi;
                                    if (x3 == x)
                                    {
                                        k_mer_bj = 1;
                                        continue;
                                    }
                                    if (x3 > x)
                                    {
                                        x3 = x;
                                    }
                                    dq_minimizer.push_back(make_pair(x3, suoyin));
                                    if (hash64(x3, mask) <= minimizer_minihash) // Same as min select new
                                    {
                                        minimizer_mini.first = x3;
                                        minimizer_mini.second = suoyin;
                                        minimizer_minihash = hash64(x3, mask);
                                        dq_min_bj = dq_minimizer.size();
                                    }
                                    if (dq_minimizer.size() < k_mer_win) // deque is incomplete
                                    {
                                        k_mer_bj = 1;
                                    }
                                }
                            }
                            else
                            {
                                x = 0;
                                x0 = 0;
                                k_mer_bj = k_mer_length;
                                dq_minimizer.clear();
                                dq_min_bj = 1;
                                minimizer_mini.first = numeric_limits<unsigned long long int>::max();
                                minimizer_mini.second = numeric_limits<unsigned long long int>::max();
                                minimizer_minihash = numeric_limits<unsigned long long int>::max();
                            }
                        }
                        else
                        {
                            continue;
                        }
                    }
                    // What are the conditions for exiting a while loop? Let's handle them
                    if (ch_atcg == '>') // Read the new sequence and start over
                    {
                        continue;
                    }
                    if (k_mer_bj == 0) // Normal exit
                    {
                        save_k_mer.push_back(make_pair(hash64(minimizer_mini.first, mask), minimizer_mini.second));
                    }
                    else // I have read to the end of the file
                    {
                        break;
                    }
                }
            }
        }
        else // If it is a line >, find a complete kmer and deque. If > is read during this process, reset it
        {
            dq_minimizer.clear();
            minimizer_mini.first = numeric_limits<unsigned long long int>::max();
            minimizer_mini.second = numeric_limits<unsigned long long int>::max();
            minimizer_minihash = numeric_limits<unsigned long long int>::max();
            dq_min_bj = 1;
            char ch_atcg;
            x = 0;
            x0 = 0;
            k_mer_bj = k_mer_length;
            while (k_mer_bj != 0 && file >> ch_atcg)
            {
                if (ch_atcg == '>')
                {
                    file.putback(ch_atcg); // Roll back one character
                    break;
                }
                if (ch_atcg != '\n' && ch_atcg != '\r')
                {
                    suoyin++;
                    if (ATUCG[ch_atcg] != 4)
                    {
                        x0 = (x0 >> 2) | ATUCG0[ch_atcg];
                        x = (x << 2) | ATUCG[ch_atcg];
                        x = (x << weiyi) >> weiyi;
                        k_mer_bj--;
                        if (k_mer_bj == 0)
                        {
                            s_kmer_sum++;
                            x3 = x0 >> weiyi;
                            if (x3 == x)
                            {
                                k_mer_bj = 1;
                                continue;
                            }
                            if (x3 > x)
                            {
                                x3 = x;
                            }
                            dq_minimizer.push_back(make_pair(x3, suoyin));
                            if (hash64(x3, mask) <= minimizer_minihash) // Same as min select new
                            {
                                minimizer_minihash = hash64(x3, mask);
                                minimizer_mini.first = x3;
                                minimizer_mini.second = suoyin;
                                dq_min_bj = dq_minimizer.size();
                            }
                            if (dq_minimizer.size() < k_mer_win) // deque is incomplete
                            {
                                k_mer_bj = 1;
                            }
                        }
                    }
                    else
                    {
                        x = 0;
                        x0 = 0;
                        k_mer_bj = k_mer_length;
                        dq_minimizer.clear();
                        dq_min_bj = 1;
                        minimizer_mini.first = numeric_limits<unsigned long long int>::max();
                        minimizer_mini.second = numeric_limits<unsigned long long int>::max();
                        minimizer_minihash = numeric_limits<unsigned long long int>::max();
                    }
                }
                else
                {
                    continue;
                }
            }
            // What are the conditions for exiting the while loop? Please go ahead and handle it.
            if (ch_atcg == '>') // Read about a new sequence and start anew.
            {
                continue;
            }
            if (k_mer_bj == 0) // Normal exit
            {
                save_k_mer.push_back(make_pair(hash64(minimizer_mini.first, mask), minimizer_mini.second));
            }
            else // Read to the end of the document
            {
                break;
            }
        }
    }
    file.close();

    unsigned long long int seed_sum = save_k_mer.size();
    // unsigned long long int max_occ = 0.0002 * s_kmer_sum + 4;
    // if (max_occ > 5000)
    //{
    // max_occ = 5000;
    //}
    // std::cout << max_occ << std::endl;

    // H-hits
    unsigned long long int seed_max = 0;

    //  save_k_mer.size();The total number of seeds
    unordered_map<unsigned long long int, unsigned long long int> map_minimizer_sum; // minimizer和其数量

    // Usage scope: for loop for iterative traversal
    for (const auto &pair : save_k_mer)
    {
        if (map_minimizer_sum.count(pair.first))
        {
            map_minimizer_sum[pair.first] = map_minimizer_sum[pair.first] + 1;
        }
        else
        {
            map_minimizer_sum.insert({pair.first, 1});
        }
    }

    for (const auto &pair : map_minimizer_sum)
    {
        seed_max = seed_max + pair.second * pair.second;
    }

    // Density and density factor
    std::cout << seed_sum / (double)s_kmer_sum << std::endl;
    std::cout << seed_sum * k_mer_win / (double)s_kmer_sum << std::endl;

    // H-hits
    std::cout << seed_max / (double)seed_sum << std::endl;

    // DR-Score
    // std::cout << seed_max / (double)seed_sum << std::endl;
    std::cout << (seed_sum * k_mer_win / (double)s_kmer_sum) * seed_max / (double)seed_sum << std::endl;
    return 0;
}
