#include<bits/stdc++.h>
#define pb push_back
#define eb emplace_back
#define ll long long
using namespace std;

struct kmer{
	int utgno;
	vector<int>rid;
};

struct contig{
	int ctind;
	string s;
	set<int>sid;
	int nxt;
	int prv;
}cntg;

int k,abd,iter2,idx,ress;

vector<contig>ctgs,ctgsn;
vector<string>edg;
vector<int>vtmp;
set<int>st;

unordered_map<string,kmer>hmap;
unordered_map<string,kmer>::iterator it;
set<int>::iterator itt;

ofstream opf("unitigs.fasta");
ofstream opf2("unitig_info.csv");
ofstream opf3("contig.fasta");

string cmap[4]={"A","T","G","C"};
string orig,tmp;

bool cmp(contig a,contig b){
	if(a.s.length()==b.s.length())
		return a.sid.size()>b.sid.size();
	return a.s.length()>b.s.length();
}

void generateUnitigs(unordered_map<string,kmer>::iterator it){
	orig=(*it).first;
	cntg.s=orig;
	cntg.sid.clear();
	idx=ctgs.size();
	tmp=orig;
	while(true){
		hmap[tmp].utgno=idx;
		edg.clear();
		for(int i=0;i<hmap[tmp].rid.size();i++){
			cntg.sid.insert(hmap[tmp].rid[i]);
		}
		hmap[tmp].rid.clear();
		for(int i=0;i<4;i++){
			if(hmap.find(tmp.substr(1,k)+cmap[i])!=hmap.end()){ 
				edg.pb(cmap[i]);
			}
		}
		if(edg.size()!=1)
			break;
		else{
			tmp=tmp.substr(1,k)+edg[0];
			if(hmap[tmp].rid.size()<1)
				break;
			for(int i=0;i<4;i++){
				if(hmap.find(cmap[i]+tmp.substr(0,k-1))!=hmap.end()){ 
					edg.pb(cmap[i]);
				}
			}
			if(edg.size()==2){
				cntg.s+=tmp[k-1];
			}
			else{
				break;
			}
		}
	}
	tmp=orig;
	while(true){
		hmap[tmp].utgno=idx;
		edg.clear();
		for(int i=0;i<hmap[tmp].rid.size();i++){
			cntg.sid.insert(hmap[tmp].rid[i]);
		}
		hmap[tmp].rid.clear();
		for(int i=0;i<4;i++){
			if(hmap.find(cmap[i]+tmp.substr(0,k-1))!=hmap.end()){ 
				edg.pb(cmap[i]);
			}
		}
		if(edg.size()!=1){
			break;
		}
		else{
			tmp=edg[0]+tmp.substr(0,k-1);
			if(hmap[tmp].rid.size()<1)
				break;
			for(int i=0;i<4;i++){
				if(hmap.find(tmp.substr(1,k)+cmap[i])!=hmap.end()){ 
					edg.pb(cmap[i]);
				}
			}
			if(edg.size()==2){
				cntg.s=edg[0]+cntg.s;
			}
			else{
				break;
			}
		}
	}
	iter2++;
	if(iter2%5000==0)
		cout<<iter2<<"\n";
	opf<<">contig "<<iter2<<"\n";
	opf<<cntg.s<<"\n";
	cntg.nxt=0;
	cntg.prv=0;
	cntg.ctind=ctgs.size();
	ctgs.eb(cntg);
	return ;
}

void make_csv(){
	opf2<<"unitig,readid_of_unitig,right_end_A,readid,right_end_T,readid,right_end_G,readid,right_end_C,readid,left_end_A,readid,left_end_T,readid,left_end_G,readid,left_end_C,readid\n";
	int idx;
	for(int i=1;i<ctgs.size();i++){
		vtmp.clear();
		int cntr=0,cntl=0;
		string strr=ctgs[i].s.substr(ctgs[i].s.length()-k+1,k-1);
		string strl=ctgs[i].s.substr(0,k-1);
		for(int j=0;j<4;j++){
			if((hmap.find(strr+cmap[j])!=hmap.end()) && (hmap[strr+cmap[j]].utgno!=i) && (hmap[strr+cmap[j]].utgno>0)){
				vtmp.eb(hmap[strr+cmap[j]].utgno);
				cntr++;
			}
			else{
				vtmp.eb(-1);
			}
		}
		for(int j=0;j<4;j++){
			if((hmap.find(cmap[j]+strl)!=hmap.end()) && (hmap[cmap[j]+strl].utgno!=i) && (hmap[cmap[j]+strl].utgno>0)){
				vtmp.eb(hmap[cmap[j]+strl].utgno);
				cntl++;
			}
			else{
				vtmp.eb(-1);
			}
		}
		if(cntr>1 || cntl>1){
			opf2<<ctgs[i].s<<",";
			for(itt=ctgs[i].sid.begin();itt!=ctgs[i].sid.end();itt++)
				opf2<<(*itt)<<"#";
			opf2<<",";
			for(int j=0;j<8;j++){
				if(vtmp[j]==-1)
					opf2<<"$,$,";
				else{
					opf2<<ctgs[vtmp[j]].s<<",";
					for(itt=ctgs[vtmp[j]].sid.begin();itt!=ctgs[vtmp[j]].sid.end();itt++)
						opf2<<(*itt)<<"#";
					opf2<<",";
				}
			}
			opf2<<"\n";
		}		
	}
}

void remove_bubbles(){
	for(int i=1;i<ctgs.size();i++)
		ctgsn.eb(ctgs[i]);
	sort(ctgsn.begin(),ctgsn.end(),cmp);
	int mxx,mxind,iter3=0;
	string suff;
	for(int i=0;i<ctgsn.size();i++){
		int ind=ctgsn[i].ctind;
		if(ctgs[ind].prv!=0 || ind==0)continue;	
		while(ctgs[ind].nxt==0){
			mxx=0;mxind=0;
			suff=ctgs[ind].s.substr(ctgs[ind].s.length()-k+1,k-1);
			for(int j=0;j<4;j++){
				if(hmap.find(suff+cmap[j])!=hmap.end()){
					if(ctgs[hmap[suff+cmap[j]].utgno].prv==0 && ctgs[hmap[suff+cmap[j]].utgno].s.length()>mxx){
						mxx=ctgs[hmap[suff+cmap[j]].utgno].s.length();
						mxind=hmap[suff+cmap[j]].utgno;
					}
				}
			}
			if(mxind==0)break;
			ctgs[ind].nxt=mxind;
			ctgs[mxind].prv=ind;
			ind=mxind;
		}
	}
	for(int i=0;i<ctgsn.size();i++){
		int ind=ctgsn[i].ctind;
		if(ctgs[ind].prv || ind==0)continue;
		opf3<<">contig "<<++iter3<<"\n";
		bool flg=false;
		while(true){
			if(flg==false){
				opf3<<ctgs[ind].s;
				flg=true;
			}
			else{
				opf3<<ctgs[ind].s.substr(k-1,5000);
			}
			ind=ctgs[ind].nxt;
			if(ind==0)
				break;
		}
		opf3<<"\n";
	}		
}

int main(int argc, char** argv){
	ios_base::sync_with_stdio(false);
	auto start=chrono::system_clock::now();
	int iter=0;
	k=stoi(argv[2]);
	abd=stoi(argv[3]);
	freopen(argv[1],"r",stdin);
	
	string read,str;
	while(getline(cin,read)){
		getline(cin,read);
		for(int i=0;i<read.length();i++){
			if(read[i]!='A' && read[i]!='T' && read[i]!='G' && read[i]!='C')
				read[i]='C';
		}
		for(int i=0;i+k<=read.length();i++){
			str=read.substr(i,k);
			hmap[str].rid.emplace_back(iter);
				
		}
		getline(cin,read);
		getline(cin,read);
		if(iter%100000==0)
			cout<<iter<<"\n";
		iter++;
	}
	cntg.s="$";
	ctgs.eb(cntg);
	cout<<"Stage 1 done\n";
	vector<unordered_map<string,kmer>::iterator>tmpp;
	cout<<"ddd "<<hmap.size()<<endl;
	for(it=hmap.begin();it!=hmap.end();it++){
		if((*it).second.rid.size()<abd){
			tmpp.eb(it);
		}
	}
	for(int i=0;i<tmpp.size();i++)
		hmap.erase(tmpp[i]);
	tmpp.clear();
	cout<<"ddd "<<hmap.size()<<endl;
	for(it=hmap.begin();it!=hmap.end();it++){
		if((*it).second.rid.size())
			generateUnitigs(it);
	}
	cout<<"Stage 2 done\n";
	//make_csv();
	remove_bubbles();
	cout<<"Stage 3 done\n";
	cout<<"ddd "<<hmap.size()<<endl;
	auto end=chrono::system_clock::now();
	chrono::duration<double>el=end-start;
	cout<<el.count();
	return 0;
}
