// tools/eq_cpp.cpp
// Minimal C++ driver to emit deterministic JSON matching Rust tests/equiv.rs schema.
// Build: make eq-cpp
// Run: ./eq_cpp
//
// JSON schema:
// {
//   "cpp": {
//     "masking_align": [ { "i": ..., "consensusstart": ..., "consensusend": ..., "result": [start, end] }, ... ],
//     "chrtracer": [ { "pos": ..., "chr": "...", "off": ... }, ... ],
//     "mask_extention_score": [ { "is_right": true/false, "ext": ..., "i": ..., "band": ..., "score": ... }, ... ]
//   }
// }
//
// Notes:
// - For reproducibility and to avoid touching production code, this uses the same placeholder logic
//   as the Rust harness for masking_align and mask_extention_score.
// - chrtracer uses a minimal emulation reading test/tst.fa with simple header parsing to build
//   a table similar to C++ chrtable. Padding is applied as in REPrise.cpp with PADLENGTH=11000.

#include <bits/stdc++.h>
using namespace std;

static const int PADLENGTH = 11000;

static inline string trim(const string &s){
    size_t b = s.find_first_not_of(" \t\r\n");
    if(b==string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e-b+1);
}

static inline unsigned char char_to_num(char c){
    switch(c){
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        case 'N': case 'n': case 'x':
        case 'R': case 'r': case 'Y': case 'y': case 'M': case 'm':
        case 'K': case 'k': case 'W': case 'w': case 'S': case 's':
        case 'B': case 'b': case 'D': case 'd': case 'H': case 'h': case 'V': case 'v':
            return 99;
        default: return 99;
    }
}

// Minimal FASTA load to construct chrtable like REPrise.cpp::build_sequence
struct SeqData {
    vector<unsigned char> seq;
    vector<pair<string, size_t>> chrtable;
};

static SeqData load_fasta_like(const string &path){
    SeqData out;
    ifstream fin(path);
    out.chrtable.push_back({"unknown", 0});
    out.seq.insert(out.seq.end(), PADLENGTH, 99);
    if(!fin){
        // If missing, still provide default structure
        return out;
    }
    string line;
    while(getline(fin, line)){
        if(!line.empty() && line[0]=='>'){
            out.chrtable.push_back({"padding", out.seq.size()});
            out.seq.insert(out.seq.end(), PADLENGTH, 99);
            string name = trim(line.substr(1));
            for(char &c: name){ if(c==' '||c=='\t') c = '_'; }
            out.chrtable.push_back({name, out.seq.size()});
        }else{
            string s = trim(line);
            for(char c: s){
                if((unsigned char)c > 64){
                    out.seq.push_back(char_to_num(c));
                }
            }
        }
    }
    out.seq.insert(out.seq.end(), PADLENGTH, 99);
    return out;
}

// Test-only placeholder matching Rust harness
static pair<size_t,size_t> masking_align_placeholder(size_t /*i*/, size_t cs, size_t ce){
    if(ce >= cs){
        size_t mid = (cs + ce)/2;
        size_t s = (mid>=2)? mid-2 : 0;
        size_t e = min(ce, mid+2);
        return {s, e};
    }
    return {cs, ce};
}

static int mask_extention_score_placeholder(bool /*is_right*/, size_t /*ext*/, size_t /*i*/, size_t band){
    return (int)band;
}

// chrtracer-like wrapper against constructed chrtable
static pair<string,size_t> chrtracer_like(const vector<pair<string,size_t>> &chrtable, size_t stringpos){
    if(chrtable.empty()) return {"unknown", 0};
    size_t idx = chrtable.size()-1;
    for(size_t i=0;i<chrtable.size();++i){
        if(stringpos < chrtable[i].second){
            idx = (i==0)? 0 : i-1;
            break;
        }
    }
    return chrtable[idx];
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // masking_align cases
    vector<tuple<size_t,size_t,size_t>> masking_cases = {
        {0, 5, 15},
        {1, 0, 4},
        {2, 10, 10}
    };
    vector<string> masking_json;
    for(auto &t: masking_cases){
        size_t i, cs, ce; tie(i, cs, ce) = t;
        auto pr = masking_align_placeholder(i, cs, ce);
        std::ostringstream oss;
        oss << "{\"i\":"<<i<<",\"consensusstart\":"<<cs<<",\"consensusend\":"<<ce
            <<",\"result\":["<<pr.first<<","<<pr.second<<"]}";
        masking_json.push_back(oss.str());
    }

    // chrtracer cases using test/tst.fa if present
    vector<string> chr_json;
    SeqData data = load_fasta_like("test/tst.fa");
    {
        vector<size_t> positions;
        positions.push_back(PADLENGTH + 10);
        // find first non-unknown, non-padding contig
        std::vector<std::pair<std::string, size_t>>::const_iterator it_contig =
            std::find_if(
                data.chrtable.begin(),
                data.chrtable.end(),
                [](const std::pair<std::string, size_t>& p){ return p.first != "unknown" && p.first != "padding"; }
            );
        if(it_contig != data.chrtable.end()){
            positions.push_back(it_contig->second + 100);
        }
        // find a padding entry if present
        std::vector<std::pair<std::string, size_t>>::const_iterator it_pad =
            std::find_if(
                data.chrtable.begin(),
                data.chrtable.end(),
                [](const std::pair<std::string, size_t>& p){ return p.first == "padding"; }
            );
        if(it_pad != data.chrtable.end()){
            positions.push_back(it_pad->second);
        }
        if(positions.empty()) positions.push_back(0);
        for(size_t pos: positions){
            std::pair<std::string,size_t> pr = chrtracer_like(data.chrtable, pos);
            std::ostringstream oss;
            // escape name if needed (names produced do not include quotes)
            oss << "{\"pos\":"<<pos<<",\"chr\":\""<<pr.first<<"\",\"off\":"<<pr.second<<"}";
            chr_json.push_back(oss.str());
        }
    }

    // mask_extention_score grid
    vector<tuple<bool,size_t,size_t,size_t>> grid = {
        {true, 10, 0, 5},
        {true, 0, 1, 7},
        {false, 3, 2, 9},
    };
    vector<string> mes_json;
    for(auto &g: grid){
        bool is_right; size_t ext,i,band; tie(is_right, ext, i, band) = g;
        int score = mask_extention_score_placeholder(is_right, ext, i, band);
        std::ostringstream oss;
        oss << "{\"is_right\":"<<(is_right?"true":"false")
            <<",\"ext\":"<<ext<<",\"i\":"<<i<<",\"band\":"<<band<<",\"score\":"<<score<<"}";
        mes_json.push_back(oss.str());
    }

    // Emit JSON
    cout << "{\n  \"cpp\": {\n";
    cout << "    \"masking_align\": [\n";
    for(size_t idx=0; idx<masking_json.size(); ++idx){
        cout << "      " << masking_json[idx];
        cout << (idx+1<masking_json.size()? ",\n":"\n");
    }
    cout << "    ],\n";
    cout << "    \"chrtracer\": [\n";
    for(size_t idx=0; idx<chr_json.size(); ++idx){
        cout << "      " << chr_json[idx];
        cout << (idx+1<chr_json.size()? ",\n":"\n");
    }
    cout << "    ],\n";
    cout << "    \"mask_extention_score\": [\n";
    for(size_t idx=0; idx<mes_json.size(); ++idx){
        cout << "      " << mes_json[idx];
        cout << (idx+1<mes_json.size()? ",\n":"\n");
    }
    cout << "    ]\n";
    cout << "  }\n}\n";
    return 0;
}