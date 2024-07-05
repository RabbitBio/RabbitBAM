//
// Created by 赵展 on 2021/4/7.
//

#include "Duplicate.h"
Duplicate::Duplicate() {
    //TODO KeyBase 严重影响性能，增大到12，会导致时间增加一倍
    this->KeyBase=12;
    this->KeyBit=1<<(2*(KeyBase));
    this->CountNum=KeyBit+100;
    this->Counts=new int[CountNum];
    memset(this->Counts,0x0,CountNum*sizeof(int));
    this->Dup=new unsigned long long[CountNum];
    this->GC=new int[CountNum];
}

Duplicate::~Duplicate() {
    delete Counts;
    delete Dup;
    delete GC;
}
uint64 Duplicate::seq2int(uint8 *b,int start,int keybit,bool &valid){
    uint64 ret = 0;
    for (int i = 0; i < keybit; i++) {
        ret<<=2;
        if ((DupvalAGCT[StatusBaseRever[bam_seqi(b,i)]&0x07]) == -1) {
            valid = false;
            return 0;
        }
        ret+=DupvalAGCT[StatusBaseRever[bam_seqi(b,i)]&0x07];
    }
    return ret;
}
void Duplicate::addRecord(uint32 key,unsigned long long kmer32,int gc){
    if (Counts[key]==0){
        //printf("ok 1\n");
        Counts[key]=1;
        Dup[key]=kmer32;
        GC[key]=gc;
    }else {
        if (Dup[key] == kmer32) {
            Counts[key]++;
            //printf("ok 2");
            //add this
            //TODO check it is still logic correct or not
            if (GC[key] > gc) GC[key] = gc;
        } else if (Dup[key] > kmer32) {
            Dup[key] = kmer32;
            Counts[key] = 1;
            GC[key] = gc;
        }
    }
}
void Duplicate::statusSeq(bam1_t *b){
    if (b->core.flag&2048) return;
    if (b->core.l_qseq<32)  return;
    int start1 = 0;
    int start2 = max(0, b->core.l_qseq - 32 - 5);

    uint8 * data = bam_get_seq(b);
    bool valid = true;

    uint64 ret = seq2int(data, start1, this->KeyBase, valid);
    uint32 key = (uint32)ret;
    //printf("key is %ud\n",key);
    if(!valid) return;

    uint64 kmer32 = seq2int(data, start2, 32, valid);
    if(!valid) return;

    int gc = 0;
    // not calculated
    if(Counts[key] == 0) {
        for(int i=0; i<b->core.l_qseq; i++) {
            if(StatusBaseRever[bam_seqi(data,i)] == 'C' || StatusBaseRever[bam_seqi(data,i)] == 'G')
                gc++;
        }
    }

    gc = round(255.0 * (double) gc / (double) b->core.l_qseq);

    addRecord(key, kmer32, (uint8)gc);
}

void Duplicate::add(Duplicate* b){
    for (int i=0;i<KeyBit;i++){
        if (b->Counts[i]&&Counts[i]){
            if (Dup[i]!=b->Dup[i]){
                if (Dup[i]>b->Dup[i]){
                    Counts[i]=b->Counts[i];
                    Dup[i]=b->Dup[i];
                    GC[i]=b->GC[i];
                }
            }else{
                Counts[i]+=b->Counts[i];
                //TODO GC的add操作不会，需要找点论文研究一下
                GC[i]=(GC[i]+b->GC[i])/2.0;
            }
        }
    }
}
double Duplicate::statAll(int* hist, double* meanGC, int histSize) {
    // histSize=32; // from Rabbit Qc
    long totalNum = 0;
    long dupNum = 0;
    int* gcStatNum = new int[histSize];
    memset(gcStatNum, 0, sizeof(int)*histSize);
    for(int key=0; key<KeyBit; key++) {
        int count = Counts[key];
        double gc = GC[key];

        if(count > 0) {
            totalNum += count;
            dupNum += count - 1;

            if(count >= histSize){
                hist[histSize-1]++;
                meanGC[histSize-1] += gc;
                gcStatNum[histSize-1]++;
            }
            else{
                hist[count]++;
                meanGC[count] += gc;
                gcStatNum[count]++;
            }
        }
    }

    for(int i=0; i<histSize; i++) {
        if(gcStatNum[i] > 0) {
            meanGC[i] = meanGC[i] / 255.0 / gcStatNum[i];
        }
    }
    delete[] gcStatNum;

    if(totalNum == 0)
        return 0.0;
    else
        return (double)(dupNum+0.0)/ totalNum;
}

