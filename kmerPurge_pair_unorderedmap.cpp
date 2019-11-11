#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/regex.hpp>
#include <tr1/unordered_map>
#include <pthread.h>

#define MAX_STRING 128
#define MODE_UNORDEREDMAP 0
#define MODE_TOKYOCABINET 1

using namespace std;
using namespace tr1;

class whiteList{
private:
	int MODE;

	ifstream *ifs1;
	ifstream *ifs2;
	ifstream *ifs3;
	ofstream *ofs1;
	ofstream *ofs2;
	unordered_map<string, int> count_map;

	int data_count;
	int min_count;
	int kmer_len;
	bool exitFlg;

	/* thread */	
	pthread_mutex_t mp1;
	pthread_mutexattr_t mattr1;
	pthread_mutex_t mp2;
	pthread_mutexattr_t mattr2;

public:
	whiteList(ifstream *ifs_a, int m, int kl){
		ifs1 = ifs_a;
		pthread_mutex_init(&mp1,NULL);
		pthread_mutex_init(&mp2,NULL);
		exitFlg = false;

		data_count = 0;
		min_count = m;
		kmer_len = kl;

		addList();
	};

	~whiteList(){
	}

	void addList(){
		string line = "";
		int freq = 0;

		while(getline(*ifs1, line)){
			if(line[0] == '>'){
				data_count++;
				if(data_count % 100000 == 0){
					printf (" <- %'d ",data_count);
					printf ("                       \r");
				}
				freq = atoi(line.substr(1).c_str());
			}else{
				if(freq >= min_count){
					count_map[line]=1;
				}
				freq = 0;
				line = "";
			}
		}
	}


	static string revComp(string line){

		reverse(line.begin(), line.end());
		for(int i=0,n=line.size();i<n;i++){
			if(line[i] == 'A'){
				line[i] = 'T';
			}else if(line[i] == 'T'){
				line[i] = 'A';
			}else if(line[i] == 'G'){
				line[i] = 'C';
			}else if(line[i] == 'C'){
				line[i] = 'G';
			}else{
				/* ATGC以外の塩基が見つかったら、エラーを出力して終了 */
				//cerr << "Nucleotide character Exception!! \"" << line[i] << "\"" << endl;
				//exit(1);
			}
		}
		return line;
	}

	static void *thread_func(void *param){
		string line;
		whiteList* s = static_cast<whiteList*>(param);

		int line_count1 = 0;
		string name1;
		string seq1;
		string third1;
		string qual1;

		int line_count2 = 0;
		string name2;
		string seq2;
		string third2;
		string qual2;

		bool outputFlg = true;

		while(!s->exitFlg){
			/* file in */		
			pthread_mutex_lock(&s->mp1);
			while(getline(*s->ifs2, line)){
				line_count1++; // 行カウントをインクリメント
				if(line_count1 == 1){ // 1行目だったら、
					name1 = line; // ヘッダを取得。
				}else if(line_count1 == 2){ // 2行目だったら、
					seq1 = line; // 塩基配列を取得。
				}else if(line_count1 == 3){ // 3行目だったら、
					third1 = line; // 3行目を取得。
				}else if(line_count1 == 4){ // 4行目だったら、
					qual1 = line; // クオリティーを取得して、
					break;
				}
			}
			while(getline(*s->ifs3, line)){
				line_count2++; // 行カウントをインクリメント
				if(line_count2 == 1){ // 1行目だったら、
					name2 = line; // ヘッダを取得。
				}else if(line_count2 == 2){ // 2行目だったら、
					seq2 = line; // 塩基配列を取得。
				}else if(line_count2 == 3){ // 3行目だったら、
					third2 = line; // 3行目を取得。
				}else if(line_count2 == 4){ // 4行目だったら、
					qual2 = line; // クオリティーを取得して、
					break;
				}
			}
			pthread_mutex_unlock(&s->mp1);


			if(line_count1!=0 && line_count2!=0){
				for(int i=0,n=(seq1.size() - s->kmer_len +1);i<n;i++){
					string f_seq = seq1.substr(i,s->kmer_len);
					string r_seq = revComp(f_seq);

					if( s->count_map.end() == s->count_map.find(f_seq) && s->count_map.end() == s->count_map.find(r_seq)){
						outputFlg = false;
						break;
					}
				}
				for(int i=0,n=(seq2.size() - s->kmer_len +1);i<n;i++){
					string f_seq = seq2.substr(i,s->kmer_len);
					string r_seq = revComp(f_seq);
					if( s->count_map.end() == s->count_map.find(f_seq) && s->count_map.end() == s->count_map.find(r_seq)){
						outputFlg = false;
						break;
					}
				}
			
				/* file out*/
				if(outputFlg){
					pthread_mutex_lock(&s->mp2);

					*s->ofs1 << name1 << endl;
					*s->ofs1 << seq1 << endl;
					*s->ofs1 << third1 << endl;
					*s->ofs1 << qual1 << endl;

					*s->ofs2 << name2 << endl;
					*s->ofs2 << seq2 << endl;
					*s->ofs2 << third2 << endl;
					*s->ofs2 << qual2 << endl;

					pthread_mutex_unlock(&s->mp2);
				}
				
				outputFlg=true;
				line_count1=0;
				line_count2=0;
			}else{
				break;
			}
		}

	}

	void makeThread(ifstream *ifs_a, ifstream *ifs_b, ofstream *ofs_a, ofstream *ofs_b, int threadNum){
		ifs2 = ifs_a;
		ifs3 = ifs_b;
		ofs1 = ofs_a;
		ofs2 = ofs_b;
		pthread_attr_t attr;
		pthread_attr_init(&attr);

		pthread_t threadN[threadNum];	

		for(int i=0 ;i<threadNum;i++){
	               	if(pthread_create(&threadN[i], NULL, thread_func, this)!=0){
               	        	perror("pthread_create()");
               		}	
		}
        	for(int i=0;i<threadNum;i++){
        	        pthread_join(threadN[i],NULL);
        	}
	}

};

string time(){ // 時間を取得する関数
	struct tm *date;
	time_t now;
	int year, month, day;
	int hour, minute, second;
	time(&now);
	date = localtime(&now);
	year = date->tm_year + 1900;
	month = date->tm_mon + 1;
	day = date->tm_mday;
	hour = date->tm_hour;
	minute = date->tm_min;
	second = date->tm_sec;

	stringstream year_str;
	year_str.precision(4);
	year_str << year;
	string year2 = year_str.str();

	stringstream month_str;
	month_str.precision(2);
	month_str << month;
	string month2 = month_str.str();

	stringstream day_str;
	day_str.precision(2);
	day_str << day;
	string day2 = day_str.str();

	stringstream hour_str;
	hour_str.precision(2);
	hour_str << hour;
	string hour2 = hour_str.str();

	stringstream minute_str;
	minute_str.precision(2);
	minute_str << minute;
	string minute2 = minute_str.str();

	stringstream second_str;
	second_str.precision(2);
	second_str << second;
	string second2 = second_str.str();

	string result = year2 + "/" + month2 + "/" + day2 + " " + hour2 + ":" + minute2 + ":" + second2;

	return result;
}

string use_popen(){ // メモリ使用量を取得する関数
	int memKB;
	int memMB;
	int memGB;
	string result;
	FILE *fp;
	char command[MAX_STRING];
	char output[MAX_STRING];
	sprintf(command, "grep VmSize /proc/%d/status", getpid());
	if ((fp = popen(command, "r")) == NULL) {
		/*Failure*/
//		return;
	}
	while (fgets(output, MAX_STRING, fp) != NULL) {
		sscanf(output,"VmSize: %d kB",&memKB);
		if(memKB > 1048576){
			memGB = memKB / 1048576;
			stringstream memGBstr;
			memGBstr << "VmSize: " << memGB << " GB";
			result = memGBstr.str();
		}else if(memKB > 1024){
			memMB = memKB / 1024;
			stringstream memMBstr;
			memMBstr << "VmSize: " << memMB << " MB";
			result = memMBstr.str();
		}else{
			stringstream memKBstr;
			memKBstr << "VmSize: " << memKB << " KB";
			result = memKBstr.str();
		}
	}
	if (pclose(fp) == -1) {
		/*Failure*/
	}
	return result;
}


int main(int argc, char *argv[]){
	setvbuf(stdout, (char *)NULL, _IONBF, 0); // 出力のバッファを無効化（進捗を上書きして出力するため）
	setlocale(LC_NUMERIC,"ja_JP.utf8"); // ロケールを設定（数字の3桁毎のカンマ区切りのため）

	string count_file = "";
	string fastq_file1 = "";
	string fastq_file2 = "";
	string out_file1   = "";
	string out_file2   = "";
	string tc_file   = "";
	int kmer_len = 17;
	int min_count = 6;
	int thread_num = 1;

	string usage = "Usage: kmerPurge -c <jf_count file> -i1 <input fastq> -i2 <input fastq> -o1 <output1 fastq> -o2 <output2 fastq> ( -k <k-mer length [17]> -m <min count [6]> -t <thread num [1]> )";

	int i;
	for(i=1; i<argc; i++){
		if(strcmp(argv[i],"-c") == 0){
			i++;
			if(argv[i] == NULL){
				cerr << usage << endl;
				exit(0);
			}else{
				count_file = argv[i];
			}
		}else if(strcmp(argv[i], "-i1") == 0){
			i++;
			if(argv[i] == NULL){
				cerr << usage << endl;
				exit(0);
			}else{
				fastq_file1 = argv[i];
			}
		}else if(strcmp(argv[i], "-i2") == 0){
			i++;
			if(argv[i] == NULL){
				cerr << usage << endl;
				exit(0);
			}else{
				fastq_file2 = argv[i];
			}
		}else if(strcmp(argv[i], "-o1") == 0){
			i++;
			if(argv[i] == NULL){
				cerr << usage << endl;
				exit(0);
			}else{
				out_file1 = argv[i];
			}
		}else if(strcmp(argv[i], "-o2") == 0){
			i++;
			if(argv[i] == NULL){
				cerr << usage << endl;
				exit(0);
			}else{
				out_file2 = argv[i];
			}
		}else if(strcmp(argv[i], "-k") == 0){
			i++;
			kmer_len = atoi(argv[i]);
		}else if(strcmp(argv[i], "-m") == 0){
			i++;
			min_count = atoi(argv[i]);
		}else if(strcmp(argv[i], "-t") == 0){
			i++;
			thread_num = atoi(argv[i]);
		}else if(strcmp(argv[i], "-help") == 0 || argc==1){
			cerr << usage << endl;
			exit(0);
		}		
	}

        if(count_file=="" || fastq_file1=="" || fastq_file2=="" || out_file1=="" || out_file2==""){
                cerr << usage << endl;
                exit(0);
        }

	cout << "Count file: " << count_file << endl;
	cout << "Fastq file1: " << fastq_file1 << endl;
	cout << "Fastq file2: " << fastq_file2 << endl;
	cout << "Output file1: " << out_file1 << endl;
	cout << "Output file2: " << out_file2 << endl;
	cout << "k-mer length: " << kmer_len << endl;
	cout << "Min of k-mer count: " << min_count << endl;
	cout << "Thread Num: " << thread_num << endl;

	/* k-mer white list open */
	ifstream ifs1(count_file.c_str());
	if(ifs1.fail()){
		cerr << "File " << count_file << " do not exist.\n";
		exit(0);
	}
	/* fastq open */
	ifstream ifs2(fastq_file1.c_str());
	if(ifs2.fail()){
		cerr << "File " << fastq_file1 << " do not exist.\n";
		exit(0);
	}
	ifstream ifs3(fastq_file2.c_str());
	if(ifs3.fail()){
		cerr << "File " << fastq_file2 << " do not exist.\n";
		exit(0);
	}
	/* out file open */
	ofstream ofs1(out_file1.c_str());
	ofstream ofs2(out_file2.c_str());

	cout << "START: reading " << count_file << " " << use_popen() << " " << time() <<  endl;
	whiteList wlist(&ifs1,min_count,kmer_len);
	cout << "END: reading " << count_file << " " << use_popen() << " " << time() <<  endl;
	
	cout << "START: reading " << fastq_file1 << " " << use_popen() << " " << time() <<  endl;
	wlist.makeThread(&ifs2,&ifs3,&ofs1,&ofs2,thread_num);
	cout << "END: reading " << fastq_file1 << " " << use_popen() << " " << time() <<  endl;
	
}

