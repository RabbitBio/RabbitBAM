//
// Created by 赵展 on 2021/4/12.
//

#include "Overrepresent.h"
Overrepresent::Overrepresent(int Total,double Center){
    this->Center=Center;
    this->Total=Total;
    this->sequence=new char*[Total+50];
    this->Pos=0;
    this->OverrepresentDate=-1.0;
}

void Overrepresent::insert(bam1_t *b) {
    if (Pos>=Total) return;
    this->sequence[Pos]=new char[MAXLEN];
    uint8* seq=bam_get_seq(b);
    this->sequence[Pos][b->core.l_qseq]='\n';
    for (int i=0;i<b->core.l_qseq;i++){
        this->sequence[Pos][i]=StatusBaseRever[bam_seqi(seq,i)];
    }
    Pos++;
}
int equals(char* x,char* y){
    for (int i=0;;i++){
        if (x[i]=='\n'&&y[i]=='\n') return 0;
        if (x[i]=='\n') return -1;
        if (y[i]=='\n') return 1;
        if (x[i]<y[i]) return -1;
        if (x[i]>y[i]) return 1;
    }
}

bool cmp(char* x,char* y){
    // return x<y
    for (int i=0;;i++){
        if (x[i]=='\n') return true;
        if (y[i]=='\n') return false;
        return x < y;
    }
}

void Overrepresent::status(){
    sort(sequence,sequence+Pos,cmp);
    printf("has been in sort\n");
    int equals_num=0;
    for (int i=0,last=0;i<Pos;i++){
        if (equals(sequence[i],sequence[last])==0){
            equals_num=max(equals_num,i-last);
        }else{
            last=i;
        }
    }
    printf("equals num is %d\n",equals_num);
    this->OverrepresentDate=(equals_num+0.0)/Pos;
}


