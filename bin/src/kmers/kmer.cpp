#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

int main( int argc, char * argv[] ) {
		if ( argc != 3 ) {
				cerr << "Usage: " << argv[0] << " fragments_size uni_fasta_sequence" << endl;
				return 1;
		}

		ifstream seqIn( argv[2] );
		map< string, size_t > peptides;
		bool first = true;
		string line, seq;
		while ( getline( seqIn, line ) ) {
				if ( line[0] == '>' ) {
						if ( !first ) {
								cerr << "   err: more than one faster sequence" << endl;
								return 1;
						}
						first = false;
						seq.clear();
				} else {
						seq += line;
				}
		}
		seqIn.close();
		if ( first ) {
				cerr << "   err: missing fasta sequence header" << endl;
				return 1;
		}
		if ( seq.empty() ) {
				cerr << "   err: empty fasta sequence" << endl;
				return 1;
		}

cout << seq << endl;

		typedef pair< size_t, vector< size_t > > countPositions_pair;
		typedef map< string, countPositions_pair > fragments_t;
		fragments_t fragments;
		size_t fragSize = atoi( argv[1] );
		for ( size_t i = 0, e = seq.size() - fragSize + 1; i < e; ++i ) {
				auto it = fragments.find( seq.substr( i, fragSize ) );
				if ( it != fragments.end() ) {
						it->second.first += 1;
						it->second.second.push_back( i );
				} else {
						fragments[seq.substr( i, fragSize )] = make_pair( 1, vector< size_t >{ i } );
				}
		}

		string filename( argv[2] );
		size_t filenameIdx = filename.find_last_of( '.' );
		string outputPrefix( filename.substr( 0, filenameIdx ) );
		ofstream fragOut( outputPrefix + ".frags" );
		for_each( begin( fragments ), end( fragments ), [&]( fragments_t::value_type const & frag ) {
				fragOut << frag.first << ' ' << frag.second.first;
				for_each( begin( frag.second.second ), end( frag.second.second ), [&]( size_t pos ) {
						fragOut << ' ' << pos;
				});
				fragOut << endl;
		});

		return 0;
}

