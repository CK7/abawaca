#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "ReadMappingReader.h"
#include "SeqIORead_fasta.h"
#include "String.h"
#include "common.h"

using namespace std;

const char* bacterial_scg = "/Users/ItaiSharon/dropbox/software/bin/scg.pl";
const char* prodigal = "/usr/local/bin/prodigal";

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Sequence {
public:
	Sequence(string _display_id, const Bio::DNAString& _seq) : display_id(_display_id), seq(_seq)	{}
	virtual 		~Sequence() 							{}
	string			get_display_id() const						{return display_id;}
	const Bio::DNAString&	get_seq() const							{return seq;}
	double			gc() const							{return Bio::gc(seq);}
	size_t			Ns() const							{return Bio::Ns(seq);}
protected:
	string			display_id;
	Bio::DNAString		seq;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Scaf_segment : public Sequence {
public:
	Scaf_segment(string _display_id, string _desc, size_t s, size_t e, const Bio::DNAString& _seq);
	string				get_desc() const					{return desc;}
	void                            add_mapped_read(const Bio::ReadMapping& rmapped, size_t dimension);
	double				get_dimension(size_t d) const				{return (d < dimensions.size())? double(dimensions[d]) : 0;}
	double				get_dimension(string dimension_name) const;
protected:
	size_t				start, end;
	string				desc;
	// Holds # of read starts mapped to this segment
	vector<double>			dimensions;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Scaf : public Sequence {
public:
	Scaf(const Bio::DNASequence& s, size_t window_size);
	~Scaf()											{for(auto it=segments.begin(); it!=segments.end(); it++) delete(*it);}
	void					add_mapped_read(const Bio::ReadMapping& rmapped, size_t dimension);
	size_t					ndps() const					{return segments.size();}
	vector<Scaf_segment*>::const_iterator	begin() const					{return segments.begin();}
	vector<Scaf_segment*>::const_iterator	end() const					{return segments.end();}
	double					cvg() const					{return double(nbps)/seq.size();}
protected:
	vector<pair<size_t, size_t> >		segment_coords;
	vector<Scaf_segment*>			segments;
	size_t					nbps = 0;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

map<string, size_t> 	dimension2index;
vector<string>		dimension_order;
size_t			this_sample;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Scaf_segment
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepare_mer_indexes_recursive(int k, string mer, map<string, size_t>& mer2index)
{
	if(k==0) {
		static int next_id = 0;
		auto it = mer2index.find(Bio::DNAString(mer).reverse_complement());
		if(it == mer2index.end()) {
			mer2index.insert(pair<string, size_t>(mer, next_id++));
			dimension_order.push_back(mer);
		}
		else {
			mer2index.insert(pair<string, size_t>(mer, it->second));
		}
		return;
	}

	static char acgt[] = {'A', 'C', 'G', 'T'};
	for(auto i=0; i<4; i++) {
		prepare_mer_indexes_recursive(k-1, mer+acgt[i], mer2index);
	}
}

void prepare_mer_indexes(int k, map<string, size_t>& mer2index)
{
	for(int i=1; i<=k; i++)
		prepare_mer_indexes_recursive(i, "", mer2index);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Scaf_segment::Scaf_segment(string _display_id, string _desc, size_t s, size_t e, const Bio::DNAString& _seq) : start(s), end(e), Sequence(_display_id, _seq), desc(_desc)
{
	int k = 4;
	if(dimension2index.size() == 0) {
		prepare_mer_indexes(k, dimension2index);
	}
	dimensions.resize(dimension_order.size());

	// Compute 1,2,3 and 4-mer frequencies
	int factor = 10;
	int N_limit = 1;
	for(auto i=1; i<=k; i++, N_limit*=10, factor*=10) {
		// Shouldn't happen...
		if(seq.size() < i)
			continue;

		unordered_map<size_t, size_t>	counts;
		int 	key=0, total=0;
		auto	it_start = seq.begin(), it_end = seq.begin();

		// Add first mer
		while(it_end-it_start < i) {
			key *= 10;	// 0=00000, 2=00010, 6=00110, 19=10011
			switch(*(it_end++)) {
				case 'A':	key += 2;	break;
				case 'C':	key += 3;	break;
				case 'G':	key += 4;	break;
				case 'T':	key += 5;	break;
				default:	key = 0;
			}
		}
		if(key > N_limit) {
			counts[key]++;
 			total++;
		}
		// Rest of the mers
		while(it_end != seq.end()) {
			key *= 10;
			key = key%factor;
			switch(*(it_end++)) {
				case 'A':	key += 2;	break;
				case 'C':	key += 3;	break;
				case 'G':	key += 4;	break;
				case 'T':	key += 5;	break;
				default:	key = 0;
			}

			// The following means that there are no N's in mer)
			if(key > N_limit) {
				counts[key]++;
				total++;
			}
		}

		for(auto it=counts.begin(); it!=counts.end(); it++) {
			string mer;
			for(int i=factor; i>0; i/=10)
			switch((it->first/i)%10) {
				case 2:	mer.push_back('A');	break;
				case 3:	mer.push_back('C');	break;
				case 4:	mer.push_back('G');	break;
				case 5:	mer.push_back('T');	break;
			}
			auto sit = dimension2index.find(mer);
			if(sit == dimension2index.end()) {
				cerr << "Error: could not find index for " << mer << endl;
				exit(-1);
			}
			dimensions[sit->second] += (total != 0)? (double(it->second)/total): 0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Scaf_segment::add_mapped_read(const Bio::ReadMapping& rmapped, size_t dimension)
{
	if(dimensions.size() <= dimension)
		dimensions.resize(dimension+1);

	size_t s = (rmapped.ref_pos() < start)? start : rmapped.ref_pos();
	size_t e = ((rmapped.ref_pos()+rmapped.read_length()-1) > end)? end : (rmapped.ref_pos()+rmapped.read_length()-1);
	dimensions[dimension] += double(e-s+1)/(this->end-this->start+1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Scaf_segment::get_dimension(string dimension_name) const
{
	auto it=dimension2index.find(dimension_name);

	return (it!=dimension2index.end())? dimensions[it->second]: 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Scaf
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Scaf::Scaf(const Bio::DNASequence& s, size_t window_size) : Sequence(s.display_id(), s.seq())
{
	// Split the scaffold into segments
	const Bio::DNAString& seq = s.seq();
	size_t			nsegments = (seq.size() - Bio::Ns(seq)) / window_size;

	if(nsegments == 0)
		nsegments = 1;
	size_t	nbps = (seq.size() - Bio::Ns(seq)) / nsegments;

	size_t	nbps_segment = 0;
	for(auto it=seq.begin(), sit=seq.begin(); it!=seq.end(); it++) {
		if(*it != 'N')
			nbps_segment++;
		if(nbps_segment == nbps) {
			Bio::DNAString subs = seq.subseq(sit-seq.begin(), it-seq.begin());
			char display_id[1024], desc[1024];

			segment_coords.push_back(pair<size_t, size_t>(sit-seq.begin()+1, it-seq.begin()+1));
			sprintf(display_id, "%s_%lu", s.display_id().c_str(), segment_coords.size());
			// SCNpilot_BF_INOC_scaffold_3677:(1,2811), 0 segment(s), 2811/2811 non-N bps
			sprintf(desc, "%s:(%lu, %lu), %lu/%lu non-Ns bps", s.display_id().c_str(), (sit-seq.begin()+1), (it-seq.begin()+1), subs.size()-Bio::Ns(subs), subs.size());
			segments.push_back(new Scaf_segment(display_id, desc, (sit-seq.begin()+1), (it-seq.begin()+1), subs));
			sit = it;
			sit++;
			nbps_segment = 0;
		}
	}

	// Last few bps may be nglected 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Scaf::add_mapped_read(const Bio::ReadMapping& rmapped, size_t dimension)
{
	size_t s=rmapped.ref_pos(), e=rmapped.ref_pos()+rmapped.read_length()-1;
	for(auto it=segment_coords.begin(); it!=segment_coords.end(); it++) {
		if(s > it->second)
			continue;
		if(e < it->first)
			break;
		segments[it-segment_coords.begin()]->add_mapped_read(rmapped, dimension);
	}

	if(dimension == this_sample)
		this->nbps += rmapped.read_length();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class Links
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Links {
public:
	struct Side {
		Side(string s, size_t e) : seq(s), end(e)				{}
		bool    operator < (const Side& other) const				{if(end==0) return false; if(other.end==0) return true; if(seq!=other.seq) {return seq<other.seq;} return end<other.end;}
		bool	operator > (const Side& other) const				{return (other < *this);}
		bool	operator ==(const Side& other) const 				{return ((seq==other.seq) && (end==other.end));}
		bool	operator !=(const Side& other) const				{return !(*this == other);}
		string	seq;
		size_t	end;
	};

	struct Link {
		Link(const Side& s1, const Side& s2) : side1(s1), side2(s2)		{if(s1>s2) {side1=s2; side2=s1;}}
		bool	operator < (const Link& other) const				{if(side1!=other.side1) return (side1<other.side1); return (side2<other.side2);}
		Side	side1, side2;
		size_t	count = 0;
	};
public:
	Links(size_t is, const map<string, Scaf*>& _scafs) : insert_size(is), scafs(_scafs)	{}
	~Links();
	void		add_pair(Bio::ReadMappingPtr mapping1, Bio::ReadMappingPtr mapping2);
	void		write_connections(ofstream& os) const;
protected:
	size_t						insert_size;
	map<Side, map<Links::Side, Links::Link*> >	links;
	const map<string, Scaf*>&			scafs;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Links::~Links()
{
	set<Links::Link*> links_to_delete;

	for(auto it=links.begin(); it!=links.end(); it++) {
		for(auto sit=it->second.begin(); sit!=it->second.end(); sit++) {
			links_to_delete.insert(sit->second);
		}
	}

	for(auto sit=links_to_delete.begin(); sit!=links_to_delete.end(); sit++)
		delete(*sit);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Links::write_connections(ofstream& os) const
{
	set<Links::Link*> links_to_write;

	for(auto it=links.begin(); it!=links.end(); it++) {
		for(auto sit=it->second.begin(); sit!=it->second.end(); sit++) {
			links_to_write.insert(sit->second);
		}
	}

	for(auto sit=links_to_write.begin(); sit!=links_to_write.end(); sit++)
		os << (*sit)->side1.seq << '\t' << (*sit)->side1.end << '\t' << (*sit)->side2.seq << '\t' << (*sit)->side2.end << '\t' << (*sit)->count << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Links::add_pair(Bio::ReadMappingPtr mapping1, Bio::ReadMappingPtr mapping2)
{
	auto scafit1 = mapping1->unmapped()? scafs.end() : scafs.find(mapping1->ref_name());
	auto scafit2 = mapping2->unmapped()? scafs.end() : scafs.find(mapping2->ref_name());

	// Either both are unmapped or both are mapped to the same sequence
	if((scafit1 == scafs.end()) && (scafit2 == scafs.end()))
		return;

	size_t end1 = 0;
	if(scafit1 != scafs.end()) {
		if((mapping1->ref_pos() < this->insert_size) && (mapping1->direction() == Bio::ReadMapping::reverse))
			end1 = 5;
		else if((mapping1->ref_pos() > scafit1->second->get_seq().size()-this->insert_size) && (mapping1->direction() == Bio::ReadMapping::forward))
			end1 = 3;
	}

	size_t end2 = 0;
	if(scafit2 != scafs.end()) {
		if((mapping2->ref_pos() < this->insert_size) && (mapping2->direction() == Bio::ReadMapping::reverse))
			end2 = 5;
		else if((mapping2->ref_pos() > scafit2->second->get_seq().size()-this->insert_size) && (mapping2->direction() == Bio::ReadMapping::forward))
			end2 = 3;
	}

	// If both reads fall into the same read this is only interesting if they connect opposite ends
	if((scafit1 == scafit2) && ((end1 == end2) || (end1 == 0) || (end2 == 0)))
		return;

	// None of the reads is mapped to the ends
	if((end1 == 0) && (end2 == 0))
		return;

	Links::Side	side1((end1 != 0)? scafit1->first : "", end1), side2((end2 != 0)? scafit2->first : "", end2);
	if((end1 != 0) && (links.find(side1) == links.end())) {
		links.insert(pair<Links::Side, map<Links::Side, Links::Link*> >(side1, map<Links::Side, Links::Link*>()));
	}
	if((end2 != 0) && (links.find(side2) == links.end())) {
		links.insert(pair<Links::Side, map<Links::Side, Links::Link*> >(side2, map<Links::Side, Links::Link*>()));
	}

	auto		it = (end1 == 0)? links.find(side2) : links.find(side1);
	auto		itm = (end1 == 0)? it->second.find(side1) : it->second.find(side2);

	if(itm == it->second.end()) {
		Links::Link* link = new Links::Link(side1, side2);
		if(end1 != 0) {
			it = links.find(side1);
			it->second.insert(pair<Links::Side, Links::Link*>(side2, link));
			itm = it->second.find(side2);
		}
		if(end2 != 0) {
			it = links.find(side2);
			it->second.insert(pair<Links::Side, Links::Link*>(side1, link));
			itm = it->second.find(side1);
		}
	}
	itm->second->count++;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main etc
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void usage(const char* prog_name)
{
	cerr << endl << "Usage: " << prog_name << " -f <fasta-file> -o <out-directory> -s <sam-files-glob> -c <source-sam-file>" << endl << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t calc_insert_size(string sam_file)
{
	Bio::SAMReader  reader(sam_file);
	size_t s = 0;
	size_t	n = 0;

	while(reader.good() && (n<20000)) {
		Bio::ReadMappingPtr mapping(reader.next_mapping());
		if(mapping == NULL) {
			break;
		}

		// Note: bowtie2 computes the wrong TLEN for inserts connecting the two ends of a circular scaffold
		if(mapping->mapped() && mapping->pair_mapped() && (mapping->insert_size() > 0) && (mapping->direction() == Bio::ReadMapping::forward)) {
			s += mapping->insert_size();
			n++;
		}
	}

	return (n>0)? s/n : 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string prepare_scgs(const string& assembly_file, const string& out_directory)
{
	char cmd[4096];

	sprintf(cmd, "mkdir %s/prodigal", out_directory.c_str());
	system(cmd);
	sprintf(cmd, "mkdir %s/scg", out_directory.c_str());
	system(cmd);

	size_t l = assembly_file.rfind('/');
	string out_name = (l == string::npos)? assembly_file : assembly_file.substr(l+1);
	char scg_name[2048], prodigal_name[2048];

	sprintf(prodigal_name, "%s/prodigal/%s", out_directory.c_str(), out_name.c_str());
	sprintf(scg_name, "%s/scg/%s", out_directory.c_str(), out_name.c_str());

	sprintf(cmd, "%s -i %s -o %s.gff -d %s.genes.fna -a %s.proteins.faa -p meta -f gff 2> %s.log", prodigal, assembly_file.c_str(), prodigal_name, prodigal_name, prodigal_name, prodigal_name);
	system(cmd);
	sprintf(cmd, "ln -s ../prodigal/%s.proteins.faa %s/scg", out_name.c_str(), out_directory.c_str());
	system(cmd);
	sprintf(cmd, "%s %s.proteins.faa 2> %s.stderr > %s.stdout", bacterial_scg, scg_name, scg_name, scg_name);
//	sprintf(cmd, "ruby %s %s.proteins.faa 2> %s.stderr > %s.stdout", bacterial_scg, scg_name, scg_name, scg_name);
	system(cmd);
	sprintf(cmd, "ln -s scg/%s.proteins.faa.bacteria.scg %s", out_name.c_str(), out_directory.c_str());
	system(cmd);

	char scg_file[2048];
	sprintf(scg_file, "%s/%s.proteins.faa.bacteria.scg", out_directory.c_str(), out_name.c_str());

	return scg_file;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char* argv[])
{
	size_t	window_size = 2000;
	size_t	max_snps = 15;

	if(argc < 5) {
		usage(argv[0]);
		return -1;
	}

	string	assembly_file;
	string	out_directory;
	string	sam_glob;
	string	this_sample_str;

	for(auto i=1; i<argc; i+=2) {
		if(!strcmp(argv[i], "-f"))
			assembly_file = argv[i+1];
		else if(!strcmp(argv[i], "-o"))
			out_directory = argv[i+1];
		else if(!strcmp(argv[i], "-s"))
			sam_glob = argv[i+1];
		else if(!strcmp(argv[i], "-c"))
			this_sample_str = argv[i+1];
		else {
			cerr << endl << "Unknown option " << argv[i] << endl << endl;
			return -1;
		}
	}

	if(assembly_file.size() == 0) {
		cerr << endl << "Error: assembly file was not specified (-f)" << endl << endl;
		return -1;
	}
	if(out_directory.size() == 0) {
		cerr << endl << "Error: output directory was not specified (-o)" << endl << endl;
		return -1;
	}
	if(sam_glob.size() == 0) {
		cerr << endl << "Error: glob path for sam files was not specified (-s)" << endl << endl;
		return -1;
	}
	if(this_sample_str.size() == 0) {
		cerr << endl << "Error: path for sam file for the assembly\'s reads was not specified (-f)" << endl << endl;
		return -1;
	}

	cerr << time_since_epoch() << "Reading assembly file (" << assembly_file << ")" << endl;
	map<string, Scaf*>			scafs;
	Bio::SeqIORead_fasta<Bio::DNASequence>  in(assembly_file);
	Bio::DNASequence*                       seq = NULL;
	size_t ndps = 0;
	while((seq = in.next_seq()) != NULL) {
		Scaf*	s = new Scaf(*seq, window_size);
 		scafs.insert(pair<string, Scaf*>(seq->display_id(), s));
		ndps += s->ndps();
	}

	char names_file[1024], lrn_file[1024], info_file[1024], links_file[1024];

	sprintf(names_file, "%s/abawaca.names", out_directory.c_str());
	sprintf(lrn_file, "%s/abawaca.lrn", out_directory.c_str());
	sprintf(info_file, "%s/abawaca.info", out_directory.c_str());
	sprintf(links_file, "%s/abawaca.links", out_directory.c_str());

	cerr << time_since_epoch() << "Creating .lrn (" << lrn_file << ") and .names (" << names_file << ") files" << endl;

	if(!directory_exists(out_directory.c_str())) {
		char	cmd[1024];
		sprintf(cmd, "mkdir %s", out_directory.c_str());
		system(cmd);
	}

	if(sam_glob.size() > 0) {
		vector<string>	sam_files = glob(sam_glob);
		for(auto vit=sam_files.begin(); vit!=sam_files.end(); vit++) {
			cerr  << time_since_epoch() << "Reading SAM file " << *vit << endl;
			dimension2index.insert(pair<string, size_t>(*vit, dimension_order.size()));
			dimension_order.push_back(*vit);

			size_t dimension = dimension2index[*vit];

			if(this_sample_str == *vit)
				this_sample = dimension;

			size_t 					insert_size = calc_insert_size(*vit);
			cerr  << time_since_epoch() << "Insert size is " << insert_size << endl;

			Bio::SAMReader  			reader(*vit);
			unordered_map<string, Bio::ReadMappingPtr>	reads;
			Links					links(insert_size, scafs);
			while(reader.good()) {
				Bio::ReadMappingPtr mapping(reader.next_mapping());
				if(mapping == NULL) {
					break;
				}

				if(dimension == this_sample) {
					string insert_name = mapping->read_name();
					size_t	p = insert_name.find_last_of('/');
					if(p != string::npos)
						insert_name.resize(p);
					auto rit = reads.find(insert_name);
					if(rit == reads.end()) {
						reads.insert(pair<string, Bio::ReadMappingPtr>(insert_name, mapping));
					}
					else {
						Bio::ReadMappingPtr mapping_other = rit->second;
						reads.erase(rit);
						links.add_pair(mapping, mapping_other);
					}
				}
				if(mapping->unmapped() || (mapping->num_snps() > max_snps) || mapping->multiple_hits()) {
					continue;
				}
				auto mit = scafs.find(mapping->ref_name());
				if(mit != scafs.end())
					mit->second->add_mapped_read(*mapping, dimension);
			}
			if(dimension == this_sample) {
				ofstream	oslinks(links_file);
				links.write_connections(oslinks);
				oslinks.close();
			}
		}
	}

	FILE* fout = fopen(names_file, "w");
	FILE* flrn = fopen(lrn_file, "w");
	FILE* finfo = fopen(info_file, "w");

	if(fout == NULL) {
		cerr << "Could not write to " << names_file << endl << endl;
		return -1;
	}
	if(flrn == NULL) {
		cerr << "Could not write to " << lrn_file << endl << endl;
		return -1;
	}
	if(finfo == NULL) {
		cerr << "Could not write to " << info_file << endl << endl;
		return -1;
	}

	fprintf(fout, "%c %lu\n", '%', ndps);
	fprintf(flrn, "%c %lu\n", '%', ndps);

	// We do not output the frequency of A/T
	fprintf(flrn, "%c %lu\n", '%', dimension_order.size());
	fprintf(flrn, "%c 9", '%');
	for(auto it=dimension_order.begin(); it!=dimension_order.end(); it++)
		if(*it != string("A"))
			fprintf(flrn, "\t1");
	fprintf(flrn, "\n");

	fprintf(flrn, "%c Key", '%');
	for(auto it=dimension_order.begin(); it!=dimension_order.end(); it++)
		if(*it != string("A"))
			fprintf(flrn, "\t%s", it->c_str());
	fprintf(flrn, "\n");

	int i=0;
	for(auto it=scafs.begin(); it!=scafs.end(); it++) {
		fprintf(finfo, "%s\t%lu\t%.3lf\t%.3lf\t%lu\n", it->second->get_display_id().c_str(), it->second->get_seq().size(), int(1000.0*it->second->cvg())/1000.0, int(1000.0*it->second->gc())/1000.0, it->second->Ns());
		for(auto sit=it->second->begin(); sit!=it->second->end(); sit++) {
			fprintf(fout, "%d\t%s\t%s\n", ++i, (*sit)->get_display_id().c_str(), (*sit)->get_desc().c_str());
			fprintf(flrn, "%d", i);
			for(auto mit=dimension_order.begin(); mit!=dimension_order.end(); mit++) {
				if(*mit != string("A"))
					fprintf(flrn, "\t%.3lf", int(1000.0*(*sit)->get_dimension(*mit))/1000.0);
			}
			fprintf(flrn, "\n");
		}
	}

	fclose(fout);
	fclose(flrn);
	fclose(finfo);

	cerr << time_since_epoch() << "Predicting genes and SCGs" << endl;
	string scg_file = prepare_scgs(assembly_file, out_directory);
	char data[2048];
	sprintf(data, "%s/data.txt", out_directory.c_str());
	FILE* fp = fopen(data, "w");

	if(fp == NULL) {
		cerr << "Could not write to " << data << endl << endl;
		return -1;
	}

	fprintf(fp, "SCG\t%s\n", scg_file.c_str());
	fprintf(fp, "Links\t%s\n", links_file);
	fprintf(fp, "Info\t%s\n", info_file);
	fprintf(fp, "Names\t%s\n", names_file);
	fprintf(fp, "Lrn\t%s\n", lrn_file);
	fprintf(fp, "Assembly\t%s\n", assembly_file.c_str());
	fclose(fp);

	cerr << time_since_epoch() << "Finished successfully" << endl;
	return 0;
}
