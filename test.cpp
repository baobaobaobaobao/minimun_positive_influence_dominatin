//����罻�������С����Ӱ��֧�伯(Positive Influence Dominating Set)���������������е�̰�Ľ����㷨������������Ƶ��µ�̰���㷨
//***Algorithm 1: ���е��㷨1(Greedy1_PIDS): һ���򵥵Ľ���̰���㷨---���Ʊ�ΪH(��),����H()�Ǻ�г����,����ͼ�������; ʱ�临�Ӷ�ΪO(n^3)
//����Ҫ����:ÿ��ѡ�������������δ����������ȵ�һ����֧����Ϊ֧���.
//��������: On positive influence dominating sets in social networks,Theoretical Computer Science 412 (2011) 265�C269
//***Algorithm 2: ���е��㷨2(Greedy2_PIDS): ��һ���򵥵�̰���㷨---û��֤�����Ʊ�, ����˵Ч��Ҫ�������̰�ķ�Ҫ��; ʱ�临�Ӷ���ΪO(n^3)
//����Ҫ����:ÿ��ѡ���������������ڵ㲻����̶Ⱥ�����һ����֧����Ϊ֧���.
//��������: A New Algorithm for Positive Influence Dominating Set in Social Networks.Proc.ASONAM 2012,253-257,IEEE

//***Algorithm 3: ��������Ƶ��㷨1(NewGreedy1_PIDS): ����̰���㷨�Ĳ��ԵĽ��(ʵ������������ܱ��ֽ������������̰�ķ����); ʱ�临�Ӷ�ΪO(n^3)
//***Algorithm 4: ��������Ƶ��㷨2(NewGreedy2_PIDS): ����̰���㷨�Ĳ��Ե���һ����ʽ�Ľ��(ʵ������������ܱ��ֽ������������̰�ķ����); ʱ�临�Ӷ�ΪO(n^3)
//***Algorithm 5: ��������Ƶ��㷨3(NewGreedy3_PIDS): �����㷨(ʵ��֤���������������̰�ķ���Ҫ��һ��,������ʱ���, �����ڴ���罻����); ʱ�临�Ӷ�ΪO(n^2)

//2016-04-16
//ע:Ϊ�˼����,�������ᵽ֧�伯DS(Dominating Set)����ָ��Ӱ��֧�伯PIDS(Positive Influence Dominating Set)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <conio.h>
#include <time.h>
#include <windows.h>
#include <math.h>

#define   GMaxVertexNum 37000             //ͼ�е������Ͻ�---����ʵ������ɵ���

//***ע:Ϊ������㷨��Ч��, ���ǶԶ�ȡ���߲�����ͼ��������Ԥ����: ���սڵ������,Ȼ�󰴶������ɴ�С���¸����нڵ���,ʹ���µ�1�Žڵ�����,�µ����ڵ����С.
//***���������㷨�����õ���Ӱ��֧�伯�ǽڵ����±�ź�Ľ��,û��ת��Ϊԭͼ�нڵ����(Ҳû�б�Ҫת����ȥ,��Ϊʵ�������ǽ���Ҫ����õ���Ӱ��֧�伯�Ĵ�С������ʱ��).


/*��־ǿ
�����Ļ�������ֻ�������Ӱ�켯��С�����������ˡ���Ϊÿ��ѡ��һ�������D�Ļ����������ÿ���������������
*/
bool      GAdjMatrix[GMaxVertexNum+1][GMaxVertexNum+1];//ͼ���ڽӾ���(0-1����),��ʾͼG.Ϊȫ�����������������ӳ�����
short int DegreeList[GMaxVertexNum+1];                 //�洢ͼ�и���Ķ�,�ʺ�����Ȳ�����2^15=32768.Ϊȫ�����������������ӳ�����
bool      DSList[GMaxVertexNum+1];                     //�洢ͼ�б�ѡ�е���Ӱ��֧�伯, ��i��Ԫ��=1��ʾ��iΪ֧���,=0��ʾ��i��֧���.Ϊȫ�����������������ӳ�����

//����1: ÿ��i�������Satisfied[i]����i������֧��������(d(i)+1)/2֮��. ֵ>=0��ʾ��iΪ�������, <0��ʾ��iΪδ������������ֵ��С��ʾδ����̶�.�ʺ�����Ȳ�����2^15=32768.Ϊȫ�����������������ӳ�����
short int Satisfied[GMaxVertexNum+1];
//����2: ÿ����֧���i(��DSList[i]=0)�ĸ�����Coverage[i]����i������δ�����ĸ���.ֵ�Ǹ�,��ʾ�õ������Ϊ֧�����ɸ�����Щ�ڵ�������.�ʺ�����Ȳ�����2^15=32768.Ϊȫ�����������������ӳ�����
short int Coverage[GMaxVertexNum+1];      //ע: һ����֧���i�����֧����(DSList[i]=1),ֵCoverage[i]��������.
//����3: ÿ����֧���i(��DSList[i]=0)�ı���Ҫ�̶�Needed[i]����i������δ�����������֮��(Ϊ��,�����ֵ��ζ�ŵ�i����Щ�ڵ���Ҫ��i��Ϊ֧����ǿ�ҳ̶�).Ϊȫ�����������������ӳ�����
int Needed[GMaxVertexNum+1];              //ע: һ����֧���i�����֧����(DSList[i]=1),ֵNeeded[i]��������. ��,Needed[..]����ЩԪ��ֵ���ܴܺ�,�����ܻ�ﵽ2*Delta*Delta, ����DeltaΪͼ������.

int    gvertexnum;                        //ͼ�е��ʵ�ʸ���
long   gedgenum;                          //ͼ�бߵ�ʵ������
short int maxdegree;                      //����ͼ�������,�ʺ�����Ȳ�����2^15=32768
short int mindegree;                      //ͼ����С��
short int evendegreenum;                  //ͼ�Ķ�Ϊż���ĵ�ĸ���
short int onedegreenum;                   //��Ϊ1�ĵ�ĸ���

//***̰�Ĳ���(�������������е�����̰���㷨��������Ƶļ����µ�̰���㷨�������õļ��ֻ�������)***

//����1: ����ѡ�񸲸����(������2: Coverage[])�ķ�֧����Ϊ֧���. ����: ����������δ�����������
//����2: ����ѡ����Ҫ�̶����(������3: Needed[])�ķ�֧����Ϊ֧���. ����: �ܻ�����Ҫ�̶����(�������)����Щδ�����ڵ�������.
//����3: ����ѡ��δ����̶����ĵ�u(������1:Satisfied[]),�ٲ����ʵ�����(�������1�����2)��ѡȡ��һ����֧������ڵ�w��Ϊ֧���, ������u��δ�����.
//����4: ����ѡ��δ����̶����ĵ�u(������1:Satisfied[]),�ٲ����ʵ�����(�������1�����2)��ѡȡ�����ɸ���֧������ڵ��Ϊ֧���, ��u��һ���Գ�Ϊ�������.


//Part0: �������
void ReadGraph(char *txtgraphfilename);    //���ļ��ж�ȡͼ�����ݣ��õ�ͼ���ڽӾ��󼰵���б�.
                                           //***ע��: Ϊ�������㷨�����Ч��,���ڸù����а��յ���ɴ�С���¸�����,ʹ�ýڵ�1�Ķ����,�ڵ�gvertexnum�Ķ���С
void CreateGraph(char *txtgraphfilename);  //�������һ��ͼ(����������ƽ����),�������ͼ��ָ�����ļ���
                                           //***ע��: Ϊ�������㷨�����Ч��, ���ڸù����а��յ���ɴ�С���¸�����,ʹ�ýڵ�1�Ķ����,�ڵ�gvertexnum�Ķ���С


//Part1: ***���������е������㷨***
//̰�Ľ����㷨1(�����з���)
int Greedy1_PIDS(void);  //̰�Ĳ���: ����1
//̰���㷨2(�����з���)
int Greedy2_PIDS(void);  //̰�Ĳ���: ����2


//Part2: ***������Ƶ����㷨***
//����Ƶ�̰���㷨1
int NewGreedy1_PIDS(void); //��������̰�Ĳ��ԵĽ��: �Ȳ��ò���1, �ڲ���1����������ѡ�����������Щ���в��ò���2����ѡ��
//����Ƶ�̰���㷨2
int NewGreedy2_PIDS(void); //��������̰�Ĳ��ԵĽ��: �Ȳ��ò���2, �ڲ���2����������ѡ�����������Щ���в��ò���1����ѡ��
//����Ƶ�̰���㷨3
int NewGreedy3_PIDS(void); //��������̰�Ĳ��ԵĽ��: �Ȳ��ò���2, �ڲ���2����������ѡ�����������Щ���в��ò���1(��Ϊ������С������)����ѡ��


//����Ƶľֲ�̰�ķ���1
int LocalGreedy1_PIDS(void);     //����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ��з���ʹ�ò���1ֱ��u����Ϊֹ
//����Ƶľֲ�̰�ķ���2
int LocalGreedy2_PIDS(void);     //����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ��з���ʹ�ò���2ֱ��u����Ϊֹ
//����Ƶľֲ�̰�ķ���3(�����㷨--���㷨ʱ�临�Ӷ�ΪO(n^2),�����㷨ʱ�临�Ӷ���O(n^3))
int LocalGreedy3_PIDS(void);     //����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ��а�����ɴ�С��������µ�֧���ֱ��u����Ϊֹ.(��ͼ��ͼʱԤ������ͼ�����Ѱ���Ƚ�������,����Ȼ���򼴿�)


//����Ƶ��µľֲ�̰�ķ���1
int NewLocalGreedy1_PIDS(void);  //����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ���ʹ�ò���1(ʹ��1��)
//����Ƶ��µľֲ�̰�ķ���2
int NewLocalGreedy2_PIDS(void);  //����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ���ʹ�ò���2(ʹ��1��)
                                 //������ķ�֧������ڵ���ʹ�ò���2(����:A����Ҫ����+BҲ�����Ҫ)
//����Ƶ��µľֲ�̰�ķ���3(�����㷨--���㷨ʱ�临�Ӷ�ΪO(n^2),�����㷨ʱ�临�Ӷ���O(n^3))
int NewLocalGreedy3_PIDS(void);  //����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ��������ȵĵ�Ϊ��֧���.(��ͼ��ͼʱԤ������ͼ�����Ѱ���Ƚ�������)



//Part4: ***�����ڸ��㷨�е��ӳ���***
int RefinePIDS(int dssize);    //�ɶ��κ��㷨�õ�����Ӱ��֧�伯DSList[1...gvertexnum]���м�С������(���СΪdssize),ȥ������Ҫ��֧������յõ���С����Ӱ��֧�伯DSList[1...gvertexnum],�������µĴ�С
bool CheckPIDS(int dssize);    //***���ڵ��Խ׶�: ���κ��㷨�õ�����Ӱ��֧�伯DSList[1...gvertexnum](���СΪdssize)������ȷ�Լ��, ��ȷ����1;���򷵻�0.

/*
//������ִ�в���(begin)*****************************************************************************************************
void main(int argc,char *argv[])
{//������ִ�иó���, ����������: ��һ�������ǿ�ִ�е��ļ���, �ڶ������������е�ͼ�ļ���, �����������Ǳ�����������ļ�

    if(argc!=3)
	{
       printf("The number of the parameters is wrong!\n");
	   return;
	}

	char *graphfilename;
	char *resultfilename;

	graphfilename=argv[1];    //��ȡ: ͼ�ļ���
	resultfilename=argv[2];   //��ȡ: ����ļ���
//������ִ�в���(end)*********************************************************************************************************

*/




//***���������е�̰�Ľ����㷨1***
int Greedy1_PIDS(void)
{//����:����ѡ�������������δ�����������(�����������)��һ����֧����Ϊ֧���.
 //ʵ�ַ���: ����ѡȡ����Coverage[...]��ȡ���ֵ���Ǹ���֧��ڵ�(DSList[]ֵ0�ĵ��Ƿ�֧���)��Ϊ֧���.
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                         //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
		//����ȳ�ʼ��Ϊ��������-1/2
		Coverage[i]=DegreeList[i];           //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)�������ڰ����е�δ����ĵ����--����¼ÿ����֧���ĸ�����
	}
    dssize=0;//֧��������ʼ��


	//Ԥ����(��): ���е��1�ĵ��ڵ�ض�ѡ��֧���(ĳЩ��Ϊ1����ҲҪ��ѡ��֧���,����һ����������1�ĵ�����ʱ)
	//̰�Ľ����㷨
	int maxcoverage, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //��ʼ��,��ʾ��ʼʱδ�����ĸ���

	while(uncoverednum)        //�������δ����Ľڵ����Ҫ��������µ�֧���
	{
		maxcoverage=0;        //��ʼ��.�ñ�����ʾ�����һ��֧�������֧����ٸ�δ����ĵ�
		for(i=1;i<=gvertexnum;i++)
			if(DSList[i]==0 && Coverage[i]>maxcoverage)  //���з�֧��㶼����Ҫѡȡ(����ĳЩ��Ϊ1�ĵ���ܱ�ѡ��֧���)
			{//̰�ĵ�ѡ����ѽڵ���Ϊ�µ�֧���
					maxcoverage=Coverage[i];
					bestvertex=i;
			}
        if(maxcoverage==0)   //��ʾ������κε���Ϊ֧��㶼���������κ�δ�����������. ����֤�����޹������ͼ������������ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //���¸ü�����Ӱ�켯�ĵ�ĸ��ڵ��δ����ȣ�ֻ��Ҫ�����Ѿ�������Ӱ�켯������ok��
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
			     if(Satisfied[j]==0)   //�����j��֧����������ܴﵽ����һ��,�����������ĵ��࣬˳����������ڵ���������㡣��õ���������ļ��ϡ�
				 {
					 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)                   //������ص�ĸ�����Coverage[.]
							if(GAdjMatrix[j][k] && DSList[k]==0)  //�ڵ��Ƿ�֧���ʱ���б�Ҫ����
								Coverage[k]--;
				 }
			}
	}


	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with Greedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

    printf("\nGreedy1_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);
    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with Greedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}




//***���ڵ��Խ׶�: ���һ����Ӱ��֧�伯����ȷ��***
bool CheckPIDS(int dssize)
{//���ڵ��Խ׶�: ���κ��㷨�õ�����Ӱ��֧�伯DSList[1...gvertexnum](���СΪdssize)������ȷ�Լ��, ��ȷ����1; ���򷵻�0.
    int i,j;
	int newdssize;
	int temp;

	//���ÿ���������Ƿ�������һ���֧���, ��ͳ��ʵ��֧���ĸ���newdssize
	newdssize=0;
	for(i=1;i<=gvertexnum;i++)
	{   if(DSList[i]==1)
	        newdssize++;
	    temp=0;
		for(j=1;j<=gvertexnum;j++)
			if(GAdjMatrix[i][j] && DSList[j]==1)
				temp++;
		if(temp<(DegreeList[i]+1)/2)  //����κε�i��������֧������������������һ��, ��ʾ�д�
			return 0;
	}
	if(newdssize!=dssize)             //���ʵ�ʵ�֧��������ԭ��������, ��ʾ�д�
		return 0;
	return 1;
}





//****��һ����Ӱ��֧�伯����һ����С��Ӱ��֧�伯****
int RefinePIDS(int dssize)
{//�ɶ��κ��㷨�õ�����Ӱ��֧�伯DSList[1...gvertexnum]���м�С������(���СΪdssize),ȥ������Ҫ��֧������յõ���С����Ӱ��֧�伯DSList[1...gvertexnum],�������µĴ�С
    int i,j,k;
	int newdssize;
	newdssize=dssize;  //��ȡԭ��Ӱ��֧�伯�Ĵ�С

	//���սڵ����С����Ĵ�����������������(ע��: ͼ�нڵ��Ѿ����ն��ɴ�С����,����1���ڵ�����,��gvertexnum�ڵ����С)
	for(i=gvertexnum;i>=1;i--)
		if(DSList[i]==1)  //��ʾDSList[i]����Ӱ�켯��һԱ���Դ˼��
		{	for(j=1;j<=gvertexnum;j++)
				if(GAdjMatrix[i][j] && Satisfied[j]==0) //��ʱ˵��֧���i�����Ը�Ϊ��֧���
					break;
			if(j==gvertexnum+1)                           //��ʱ˵��֧���i���κ��ڵ�j����Satisfied[j]<0, ��˽ڵ�i���Ը�Ϊ��֧���
			{
				newdssize--;
				DSList[i]=0;
		        for(k=1;k<=gvertexnum;k++)
			        if(GAdjMatrix[i][k])
			           Satisfied[k]--;
			}
		}

	return newdssize;//���ؼ�С�������Ӱ��֧�伯�Ĵ�С
}





/*���ڿ�ʼ����ʦ�����ġ�
��һ����ʵ��Crossver operator*/



void  CrossverOperator()
{

/*���⣬���x��y��ζ��ƣ�*/

    bool      DSList_Z[GMaxVertexNum+1];       //���Z
    bool      DSList_X[GMaxVertexNum+1];      //���X��
    bool      DSList_Y[GMaxVertexNum+1];      //���Y��
//����������x��y,һ����x��һ����y��
    for(int i=0; i<(strlen(DSList)/2); i++)
    {
        DSList_X[i]=DSList[i];
        DSList_Y[i]=DSList[(strlen(DSList)/2)+i];
    }

    //��T��Ϊ�ռ���
    for(int i=0;i<(strlen(DSList);i++){
         DSList[i]=0;
    }
    /*��x��y���н��档����z������T*/
    for(int i=0; i<strlen(DSList_X); i++)
    {

        if(DSList_X[i]==DSList_Y[i]&&DSList_X[i]==1)
        {
            DSList_Z[i]=1;

        }
        DSList_Z[i]=0;
    }

    /*��z���в�����*/

    for(int i=0; i<strlen(DSList_X); i++)
    {
        if(DSList_Z[i]==1)
        {
            //����һ��0��1���������
            float RangNum=rand()/(RAND_MAX+1.0);
            if(RangNum<0.15)
            {
                DSList_Z[i]=0;

            }
        }
    }


   printf("\n������Z������: %5d\n",strlen(DSList_Z));
}





































//***���������е�̰���㷨2***
int Greedy2_PIDS(void)
{//����: ��������������������ڵ�δ����̶Ⱥ����(������Ҫ�̶����)��һ����֧����Ϊ֧���,ֱ�����е㶼����Ϊֹ.
 //ʵ�ַ���: ����ѡȡ����Needed[...]��ȡֵ��С(������ֵ���)���Ǹ���֧��ڵ�(DSList[.]ֵ0�ĵ��Ƿ�֧���)��Ϊ֧���.
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                       //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
	}

    for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		Needed[i]=0;         //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)������δ�����ڵ��δ�����֮��(Ϊ����,�����ֵ��ʾ�÷�֧��㱻��Ҫ�ĳ̶�)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//֧��������ʼ��


	//Ԥ����(��): ���е��1�ĵ��ڵ�ض�ѡ��֧���(ĳЩ��Ϊ1����ҲҪ��ѡ��֧���,����һ����������1�ĵ�����ʱ)
	//̰���㷨
	int minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //��ʼ��,��ʾ��ʼʱδ�����ĸ���

	while(uncoverednum)        //�������δ����Ľڵ����Ҫ��������µ�֧���
	{
		minneeded=0;           //��ʼ��Ϊ0.һ��Ϊ����
		for(i=1;i<=gvertexnum;i++)
			if(DSList[i]==0 && Needed[i]<minneeded)  //���з�֧��㶼����Ҫѡȡ(����ĳЩ��Ϊ1�ĵ���ܱ�ѡ��֧���)
			{//̰�ĵ�ѡ����ѽڵ���Ϊ�µ�֧���
					minneeded=Needed[i];
					bestvertex=i;
			}
        if(minneeded==0)   //��ʾ������κε���Ϊ֧��㶼���������κ�δ�����������. ����֤�����޹������ͼ������������ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //���¸��ڵ��δ����̶�
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
  		         if(Satisfied[j]<=0)   //ֻ�ж�֮ǰ��δ������ڵ�(�ոջ������ĵ�֮ǰҲ��δ�����)���б�Ҫ���������Ϣ
				 {	 for(k=1;k<=gvertexnum;k++)               //������ص���ڵ�Ĳ�����֮��Needed[]
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //�ڵ��Ƿ�֧���ʱ���б�Ҫ����
								Needed[k]++;                  //ע: ֵΪ��,�����ֵ��ʾ����Ҫ�ĳ̶�
				 }
				 if(Satisfied[j]==0)                          //�����j��֧����������ܴﵽ����һ��,�����������ĵ���
						 uncoverednum--;
			}
	}

    //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with Greedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);
    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with Greedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}


//***����Ƶ�̰���㷨1***
int NewGreedy1_PIDS(void)
{//���㷨ʵ����������������̰���㷨���Ե�һ�ֽ�ϵ�Ӧ��.
 //����1:����ѡ������������δ�����������(�����������)��һ����֧����Ϊ֧���.
 //ʵ�ַ���: ����ѡȡCoverage[...]��ȡ���ֵ���Ǹ���֧��ڵ�u(DSList[u]=0)��Ϊ֧���.
 //����2: ����ж�����������Ľڵ�u��ʱ, ����ѡ���Ǹ����������ڵ㲻����̶Ⱥ;���ֵ���(������Ҫ�̶����)��һ����֧���u��Ϊ֧���.
 //ʵ�ַ���: ��Coverage[u]ֵ���Ķ���ڵ�u������ѡȡNeeded[u]ֵ��С(������ֵ���)���Ǹ���֧��ڵ�u(DSList[u]=0)��Ϊ֧���.
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                         //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
		Coverage[i]=DegreeList[i];           //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)�������ڰ����е�δ����ĵ����--����¼ÿ����֧���ĸ�����
	}

    for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		Needed[i]=0;                 //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)������δ�����ڵ��δ�����֮��(Ϊ����,�����ֵ��ʾ�÷�֧��㱻��Ҫ�ĳ̶�)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//֧��������ʼ��


	//Ԥ����(��): ���е��1�ĵ��ڵ�ض�ѡ��֧���(ĳЩ��Ϊ1����ҲҪ��ѡ��֧���,����һ����������1�ĵ�����ʱ)
	//�Ľ���̰�Ľ����㷨
	int maxcoverage, minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //��ʼ��,��ʾ��ʼʱû�б������ĸ���

	while(uncoverednum)        //�������û�б�����Ľڵ����Ҫ��������µ�֧���
	{
		maxcoverage=0;         //��ʼ��.����̰�Ĳ���1
		minneeded=0;           //��ʼ��.����̰�Ĳ���2
		for(i=1;i<=gvertexnum;i++)
		{	if(DSList[i]==0)   //���з�֧��㶼����Ҫѡȡ(����ĳЩ��Ϊ1�ĵ������ѡ��֧���)
			{
				if(Coverage[i]>maxcoverage)  //̰�Ĳ���1: ����ѡ�񸲸����ķ�֧���Ϊ�µ�֧���
				{
					maxcoverage=Coverage[i];
					bestvertex=i;
				}
				else if(Coverage[i]==maxcoverage && Coverage[i]>0 && Needed[i]<minneeded)
				{//̰�Ĳ���2: ������ͬ�����������������ѡ����Ҫ�̶����(��Satisfied[u]��С��u)�ķ�֧���Ϊ�µ�֧���
					minneeded=Needed[i];
					bestvertex=i;
				}
			}
		}
        if(maxcoverage==0)   //��ʾ������κε���Ϊ֧��㶼���������κ�û������ĵ�ı�֧��̶�. ����֤���������Ӧ�ò��ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //���¸��ڵ�δ����̶�
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
				 if(Satisfied[j]<=0) //ֻ�ж�֮ǰ��δ������ڵ�(�ոջ������ĵ�֮ǰҲ��δ�����)���б�Ҫ���������Ϣ
				 {	 for(k=1;k<=gvertexnum;k++)               //������ص���ڵ�Ĳ�����֮��Needed[.]**************
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //�ڵ��Ƿ�֧���ʱ���б�Ҫ����*************
								Needed[k]++;          //ע�⣺ֵΪ���������ֵ��ʾ����Ҫ�ĳ̶�*************
				 }
			     if(Satisfied[j]==0) //�����j��֧����������ܴﵽ����һ��,���������ĵ���
				 {	 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)               //������ص�ĸ�����Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)  //*************
							Coverage[k]--;                    //*************
				 }
			}
	}

   //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewGreedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�


    printf("\nNewGreedy1_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);
	//�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewGreedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewGreedy1_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}



//***����Ƶ�̰���㷨2***
int NewGreedy2_PIDS(void)
{//���㷨ʵ����������������̰���㷨���Ե�һ�ֽ�ϵ�Ӧ��.
 //����2: ����ѡ�����������ڵ㲻����̶Ⱥ;���ֵ���(������Ҫ�̶����)��һ����֧���u��Ϊ֧���.
 //ʵ�ַ���: ����ѡȡNeeded[u]ֵ��С(������ֵ���)���Ǹ���֧��ڵ�u(DSList[u]=0)��Ϊ֧���.
 //����1: ����ж������Ҫ�̶����Ľڵ�ʱ, ������ѡ������������δ�����������(�����������)��һ����֧���u��Ϊ֧���.
 //ʵ�ַ���: ��Needed[u]����ֵ���Ķ���ڵ�u������ѡȡCoverage[...]���ֵ���Ǹ���֧��ڵ�u(DSList[u]=0)��Ϊ֧���.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                         //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
		Coverage[i]=DegreeList[i];           //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)�������ڰ����е�δ����ĵ����--����¼ÿ����֧���ĸ�����
	}

    for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		Needed[i]=0;                 //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)������δ�����ڵ��δ�����֮��(Ϊ����,�����ֵ��ʾ�÷�֧��㱻��Ҫ�ĳ̶�)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//֧��������ʼ��


	//Ԥ����(��): ���е��1�ĵ��ڵ�ض�ѡ��֧���(ĳЩ��Ϊ1����ҲҪ��ѡ��֧���,����һ����������1�ĵ�����ʱ)
	//�Ľ���̰�Ľ����㷨
	int maxcoverage, minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //��ʼ��,��ʾ��ʼʱû�б������ĸ���

	while(uncoverednum)        //�������û�б�����Ľڵ����Ҫ��������µ�֧���
	{
		maxcoverage=0;         //��ʼ��.����̰�Ĳ���1
		minneeded=0;           //��ʼ��.����̰�Ĳ���2
		for(i=1;i<=gvertexnum;i++)
		{	if(DSList[i]==0)   //���з�֧��㶼����Ҫѡȡ(����ĳЩ��Ϊ1�ĵ������ѡ��֧���)
			{
				if(minneeded>Needed[i])  //̰�Ĳ���2: ����ѡ����Ҫ�̶����(��Needed[u]��С��u)�ķ�֧���Ϊ�µ�֧���
				{
					minneeded=Needed[i];
					bestvertex=i;
				}
				else if(Needed[i]==minneeded && Needed[i]<0 && Coverage[i]>maxcoverage)
				{//̰�Ĳ���1: ����Ҫ�̶�ͬ�����������������ѡ�񸲸�����һ����֧���Ϊ�µ�֧���
					maxcoverage=Coverage[i];
					bestvertex=i;
				}
			}
		}
        if(minneeded==0)   //��ʾ������κε���Ϊ֧��㶼���������κ�û������ĵ�ı�֧��̶�. ����֤���������Ӧ�ò��ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //���¸��ڵ�δ����̶�
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
				 if(Satisfied[j]<=0)   //ֻ�ж�֮ǰ��δ������ڵ�(�ոջ������ĵ�֮ǰҲ��δ�����)���б�Ҫ���������Ϣ
				 {	 for(k=1;k<=gvertexnum;k++)               //������ص���ڵ�Ĳ�����֮��Needed[.]**************
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //�ڵ��Ƿ�֧���ʱ���б�Ҫ����*************
								Needed[k]++;                  //ע��ֵΪ���������ֵ��ʾ����Ҫ�ĳ̶�*************
				 }
			     if(Satisfied[j]==0) //�����j��֧����������ܴﵽ����һ��,���������ĵ���
				 {	 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)               //������ص�ĸ�����Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)  //*************
							Coverage[k]--;                    //*************
				 }
			}
	}

   //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewGreedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

    printf("\nNewGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);
    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewGreedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}


//***����Ƶ�̰���㷨3***
int NewGreedy3_PIDS(void)
{//���㷨ʵ����������������̰���㷨���Ե�һ�ֽ�ϵ�Ӧ��.
 //����2: ����ѡ�����������ڵ㲻����̶Ⱥ;���ֵ���(������Ҫ�̶����)��һ����֧���u��Ϊ֧���.
 //ʵ�ַ���: ����ѡȡNeeded[u]ֵ��С(������ֵ���)���Ǹ���֧��ڵ�u(DSList[u]=0)��Ϊ֧���.
 //����1: ����ж������Ҫ�̶����Ľڵ�ʱ, ������ѡ����������С��δ�����������(���������� С)��һ����֧���u��Ϊ֧���.
 //ʵ�ַ���: ��Needed[u]����ֵ���Ķ���ڵ�u������ѡȡCoverage[...]��С ֵ(��ʱƽ��ÿ��δ�����Ĳ�����Ⱥܴ�)���Ǹ���֧��ڵ�u(DSList[u]=0)��Ϊ֧���.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                         //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2;   //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
		Coverage[i]=DegreeList[i];           //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)�������ڰ����е�δ����ĵ����--����¼ÿ����֧���ĸ�����
	}

    for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		Needed[i]=0;                 //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)������δ�����ڵ��δ�����֮��(Ϊ����,�����ֵ��ʾ�÷�֧��㱻��Ҫ�ĳ̶�)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//֧��������ʼ��


	//Ԥ����(��): ���е��1�ĵ��ڵ�ض�ѡ��֧���(ĳЩ��Ϊ1����ҲҪ��ѡ��֧���,����һ����������1�ĵ�����ʱ)
	//�Ľ���̰�Ľ����㷨
	int mincoverage, minneeded, bestvertex;
	int uncoverednum;
	uncoverednum=gvertexnum;   //��ʼ��,��ʾ��ʼʱû�б������ĸ���

	while(uncoverednum)        //�������û�б�����Ľڵ����Ҫ��������µ�֧���
	{
		mincoverage=gvertexnum;  //��ʼ��.����̰�Ĳ���1
		minneeded=0;             //��ʼ��.����̰�Ĳ���2
		for(i=1;i<=gvertexnum;i++)
		{	if(DSList[i]==0)   //���з�֧��㶼����Ҫѡȡ(����ĳЩ��Ϊ1�ĵ������ѡ��֧���)
			{
				if(minneeded>Needed[i])  //̰�Ĳ���2: ����ѡ����Ҫ�̶����(��Needed[u]��С��u)�ķ�֧���Ϊ�µ�֧���
				{
					minneeded=Needed[i];
					bestvertex=i;
				}
				else if(Needed[i]==minneeded && Needed[i]<0 && Coverage[i]<mincoverage)
				{//̰�Ĳ���1: ����Ҫ�̶�ͬ�����������������ѡ�񸲸�����һ����֧���Ϊ�µ�֧���
					mincoverage=Coverage[i];
					bestvertex=i;
				}
			}
		}
        if(minneeded==0)   //��ʾ������κε���Ϊ֧��㶼���������κ�û������ĵ�ı�֧��̶�. ����֤���������Ӧ�ò��ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}
		DSList[bestvertex]=1;          //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //���¸��ڵ�δ����̶�
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
				 if(Satisfied[j]<=0)   //ֻ�ж�֮ǰ��δ������ڵ�(�ոջ������ĵ�֮ǰҲ��δ�����)���б�Ҫ���������Ϣ
				 {	 for(k=1;k<=gvertexnum;k++)               //������ص���ڵ�Ĳ�����֮��Needed[.]**************
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //�ڵ��Ƿ�֧���ʱ���б�Ҫ����*************
								Needed[k]++;                  //ע��ֵΪ���������ֵ��ʾ����Ҫ�ĳ̶�*************
				 }
			     if(Satisfied[j]==0) //�����j��֧����������ܴﵽ����һ��,���������ĵ���
				 {	 uncoverednum--;
					 for(k=1;k<=gvertexnum;k++)               //������ص�ĸ�����Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)  //*************
							Coverage[k]--;                    //*************
				 }
			}
	}

   //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewGreedy3_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

    printf("\nNewGreedy3_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);
    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewGreedy3_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewGreedy3_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}



//***����Ƶľֲ�̰�ķ���1***
int LocalGreedy1_PIDS(void)
{//����: ÿ��ѡ��һ��������u, Ȼ����u�ķ�֧������ڵ��з���ʹ�ò���1ֱ��u����Ϊֹ.
 //������u��Satisfied[.]ֵΪ���Ҿ���ֵ���ĵ�.�ڵ�u����������������ѡ�񸲸����������ɷ�֧����Ϊ֧���(������ѡ��Coverage[.]ֵ������u���ڵķ�֧���),ֱ����u����Ϊֹ.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                       //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
		Coverage[i]=DegreeList[i];         //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)�������ڰ����е�δ����ĵ����--����¼ÿ����֧���ĸ�����
	}
    dssize=0;//֧��������ʼ��

	//�ֲ�̰���㷨
	int uncoverednum;
	int currentvertex, minsatisfied, maxcoverage, maxcvertex;

	uncoverednum=gvertexnum;   //��ʼ��,��ʾ��ʼʱδ�����ĸ���
	while(uncoverednum)
	{//���uncoverednum>0, ������δ����Ľڵ����Ҫ��������µ�֧���

        minsatisfied=0;        //�����ҳ�����̶���С�Ľڵ�(��Satisfied[]ֵΪ����ֵ��С�ĵ�,�����ֵ��ʾ��δ����̶ȵ����)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //Ѱ�ҵ�ǰ����̶���С�Ľڵ㴦��(����򵥵㴦��,Ҳ���԰��������δ����ڵ����ȴ���)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******�ֲ�̰�Ĳ���:һ���Խ���ǰδ����Ľڵ�currentvertex��Ϊ����ڵ�(���ν��������и������������ɷ�֧����Ϊ֧���)
		while(Satisfied[currentvertex]<0) //Satisfied[currentvertex]<0��ʾcurrentvertex��δ�����
		{
			maxcoverage=0;   //�����ҳ��ܸ���δ����ڵ�����һ����֧���(currentvertex���ڵ�)

			//����ȡ�����漸�б��Ϊ  //**�Ĵ���(���ʵ�鷢�ֲ����ǵ�ǰ��Ҳ�ɱ�Ϊ֧���ʱЧ�����õĻ�)
			//if(DSList[currentvertex]==0 && Coverage[currentvertex]>maxcoverage)  //** �����ǰδ����ڵ�currentvertex�Ƿ�֧���ʱ,Ҳ�������ɱ�Ϊ֧���
			//{	maxcoverage=Coverage[currentvertex];  //*
			//	maxcvertex=currentvertex;      //*
			//}

			for(i=1;i<=gvertexnum;i++)
			{
				if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Coverage[i]>maxcoverage)  //�ֲ�̰�ĵ�ѡ��һ���ڵ�(��֧���)
				{   	maxcoverage=Coverage[i];
						maxcvertex=i;
				}
			}
		    if(maxcoverage==0)   //��ʾ������κ��ڵ���Ϊ֧��㶼���������κ�δ����ĵ������̶�. ����֤����������ڲ����������ͼ�в��ᷢ��
			{
				printf("No solution to the instance!\n");
				return 0;
			}

			DSList[maxcvertex]=1;          //�����µ�֧���
			dssize++;                      //֧����������
			for(j=1;j<=gvertexnum;j++)     //������ص��δ����̶�Satisfied[.]
				if(GAdjMatrix[maxcvertex][j])
				{   Satisfied[j]++;
					if(Satisfied[j]==0)    //�����j��֧����������ܴﵽ����һ��,���������ĵ���
					{	uncoverednum--;
					    for(k=1;k<=gvertexnum;k++)               //������ص�ĸ�����Coverage[.]**************
							if(GAdjMatrix[j][k] && DSList[k]==0)
								Coverage[k]--;
					}
				}
		}
	}


	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with LocalGreedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nLocalGreedy1_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);
	//�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with LocaGreedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nLocalGreedy1_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}



//***����Ƶľֲ�̰�ķ���2***
int LocalGreedy2_PIDS(void)
{//����: ÿ��ѡ��һ��������u, Ȼ����u�ķ�֧������ڵ��з���ʹ�ò���2ֱ��u����Ϊֹ.
 //������u��Satisfied[.]ֵΪ���Ҿ���ֵ���ĵ�.�ڵ�u����������������ѡ����Ҫ�̶��������ɷ�֧����Ϊ֧���(������ѡ��Needed[.]����ֵ������u���ڵķ�֧���),ֱ����u����Ϊֹ.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                       //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
	}
    for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		Needed[i]=0;         //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)������δ�����ڵ��δ�����֮��(Ϊ����,�����ֵ��ʾ�÷�֧��㱻��Ҫ�ĳ̶�)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//֧��������ʼ��

	//�ֲ�̰���㷨
	int uncoverednum;
	int currentvertex, minsatisfied, minneeded, minneededvertex;

	uncoverednum=gvertexnum;         //��ʼ��,��ʾ��ʼʱδ�����ĸ���
	while(uncoverednum)
	{//���uncoverednum>0, ������δ����Ľڵ����Ҫ��������µ�֧���

        minsatisfied=0;              //�����ҳ�����̶���С�Ľڵ�(��Satisfied[]ֵΪ����ֵ��С�ĵ�,�����ֵ��ʾ��δ����̶ȵ����)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //Ѱ�ҵ�ǰ����̶���С�Ľڵ㴦��(����򵥵㴦��,Ҳ���԰��������δ����ڵ����ȴ���)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******�ֲ�̰�ķ�����:һ���Խ���ǰδ����Ľڵ�currentvertex��Ϊ����ڵ�(���ν��������б���Ҫ�̶��������ɷ�֧����Ϊ֧���)
		while(Satisfied[currentvertex]<0) //Satisfied[currentvertex]<0��ʾcurrentvertex��δ�����
		{
			minneeded=0;   //�����ҳ�����Ҫ�̶����(��Needed[]ֵΪ���Ҿ���ֵ���)��һ����֧���(currentvertex���ڵ�)

			//����ȡ�����漸�б��Ϊ  //**�Ĵ���(���ʵ�鷢�ֲ����ǵ�ǰ��Ҳ�ɱ�Ϊ֧���ʱЧ�����õĻ�)
			//if(DSList[currentvertex]==0 && Needed[currentvertex]<minneeded)  //** �����ǰδ����ڵ�currentvertex�Ƿ�֧���ʱ,Ҳ�������ɱ�Ϊ֧���
			//{	minneeded=Needed[currentvertex];  //*
			//	minneededvertex=currentvertex;      //*
			//}

			for(i=1;i<=gvertexnum;i++)
			{
				if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Needed[i]<minneeded)  //�ֲ�̰�ĵ�ѡ��һ���ڵ�(��֧���)
				{   	minneeded=Needed[i];
						minneededvertex=i;
				}
			}
		    if(minneeded==0)   //��ʾ������κ��ڵ���Ϊ֧��㶼���������κ�δ����ĵ������̶�. ����֤����������ڲ����������ͼ�в��ᷢ��
			{
				printf("No solution to the instance!\n");
				return 0;
			}

			DSList[minneededvertex]=1;          //�����µ�֧���
			dssize++;                      //֧����������
			for(j=1;j<=gvertexnum;j++)     //������ص��δ����̶�Satisfied[.]
				if(GAdjMatrix[minneededvertex][j])
				{    Satisfied[j]++;
  					if(Satisfied[j]<=0)   //ֻ�ж�֮ǰ��δ������ڵ�(�ոջ������ĵ�֮ǰҲ��δ�����)���б�Ҫ���������Ϣ
					{	for(k=1;k<=gvertexnum;k++)                //������ص���ڵ�Ĳ�����֮��Needed[]
    						if(GAdjMatrix[j][k] && DSList[k]==0)  //�ڵ��Ƿ�֧���ʱ���б�Ҫ����
								Needed[k]++;                      //ע: ֵΪ��,�����ֵ��ʾ����Ҫ�ĳ̶�
					}
					if(Satisfied[j]==0)                           //�����j��֧����������ܴﵽ����һ��,�����������ĵ���
						 uncoverednum--;
				}
		}
	}

    //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with LocalGreedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nLocalGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);

    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with LocaGreedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nLocalGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}


//***����Ƶľֲ�̰�ķ���3(�����㷨--���㷨ʱ�临�Ӷ�ΪO(n^2),�����㷨ʱ�临�Ӷ���O(n^3))***
int LocalGreedy3_PIDS(void)
{//����: ÿ��ѡ��һ��������u, Ȼ����u�ķ�֧������ڵ��з�����ӵ�ȴ�ĵ���Ϊ�µ�֧���,ֱ��u����Ϊֹ.
 //������u��Satisfied[.]ֵΪ���Ҿ���ֵ���ĵ�.��u������������ӵ�ȴ�ĵ���Ϊ�µ�֧���(��ͼ��ͼʱԤ������ͼ�����Ѱ���Ƚ�������,����Ȼ���򼴿�),ֱ����u����Ϊֹ.

	int i,j,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                       //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
	}

    dssize=0;//֧��������ʼ��

	//�ֲ�̰���㷨
	int uncoverednum;
	int currentvertex, minsatisfied, bestvertex;

	uncoverednum=gvertexnum;         //��ʼ��,��ʾ��ʼʱδ�����ĸ���
	while(uncoverednum)
	{//���uncoverednum>0, ������δ����Ľڵ����Ҫ��������µ�֧���

        minsatisfied=0;              //�����ҳ�����̶���С�Ľڵ�(��Satisfied[]ֵΪ����ֵ��С�ĵ�,�����ֵ��ʾ��δ����̶ȵ����)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //Ѱ�ҵ�ǰ����̶���С�Ľڵ㴦��(����򵥵㴦��,Ҳ���԰��������δ����ڵ����ȴ���)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******�ֲ�̰�Ĳ���:һ���Խ���ǰδ����Ľڵ�currentvertex��Ϊ����ڵ�(���ν��������е���������ɷ�֧����Ϊ֧���)
		while(Satisfied[currentvertex]<0) //Satisfied[currentvertex]<0��ʾcurrentvertex��δ�����
		{
			for(i=1;i<=gvertexnum;i++)
				if(GAdjMatrix[currentvertex][i] && DSList[i]==0)  //�ֲ�̰�ĵ�ѡ��һ���ڵ�(��֧���)--��ȴ�����(��Ȼ���򼴱�ʾ�ȴ������,������Ϊͼ�����Ѱ���Ƚ�������)
				{		bestvertex=i;
				        break;
				}
		    if(i==gvertexnum+1)   //��ʾ������κ��ڵ���Ϊ֧��㶼���������κ�δ����ĵ������̶�. ����֤����������ڲ����������ͼ�в��ᷢ��
			{
				printf("No solution to the instance!\n");
				return 0;
			}

			DSList[bestvertex]=1;          //�����µ�֧���
			dssize++;                      //֧����������
			for(j=1;j<=gvertexnum;j++)     //������ص��δ����̶�Satisfied[.]
				if(GAdjMatrix[bestvertex][j])
				{    Satisfied[j]++;
  					 if(Satisfied[j]==0)   //�����j��֧����������ܴﵽ����һ��,�����������ĵ���
						 uncoverednum--;
				}
		}
	}


    //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with LocalGreedy3_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nLocalGreedy3_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);

    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with LocaGreedy3_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nLocalGreedy3_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}




//***����Ƶ��µľֲ�̰�ķ���1***
int NewLocalGreedy1_PIDS(void)
{//����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ���ʹ�ò���1(ʹ��1��)
 //������u��Satisfied[.]ֵΪ���Ҿ���ֵ���ĵ�.�ڵ�u��������ѡ�񸲸������(��Coverage[.]ֵ������u����)��һ����֧����Ϊ֧���.

	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                       //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
		Coverage[i]=DegreeList[i];         //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)�������ڰ����е�δ����ĵ����--����¼ÿ����֧���ĸ�����
	}
    dssize=0;//֧��������ʼ��

	//�µľֲ�̰���㷨
	int uncoverednum;
	int currentvertex, minsatisfied, maxcoverage, maxcvertex;

	uncoverednum=gvertexnum;   //��ʼ��,��ʾ��ʼʱδ�����ĸ���
	while(uncoverednum)
	{//���uncoverednum>0, ������δ����Ľڵ����Ҫ��������µ�֧���

        minsatisfied=0;        //�����ҳ�����̶���С�Ľڵ�(��Satisfied[]ֵΪ����ֵ��С�ĵ�,�����ֵ��ʾ��δ����̶ȵ����)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //Ѱ�ҵ�ǰ����̶���С�Ľڵ㴦��(����򵥵㴦��,Ҳ���԰��������δ����ڵ����ȴ���)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******�µľֲ�̰�Ĳ���: ����ǰδ����ڵ�currentvertex�������и���������һ����֧����Ϊ֧���
		maxcoverage=0;   //�����ҳ��ܸ���δ����ڵ�����һ����֧���(currentvertex���ڵ�)

		//����ȡ�����漸�б��Ϊ  //**�Ĵ���(���ʵ�鷢�ֲ����ǵ�ǰ��Ҳ�ɱ�Ϊ֧���ʱЧ�����õĻ�)
		//if(DSList[currentvertex]==0 && Coverage[currentvertex]>maxcoverage)  //** �����ǰδ����ڵ�currentvertex�Ƿ�֧���ʱ,Ҳ�������ɱ�Ϊ֧���
		//{	maxcoverage=Coverage[currentvertex];  //*
		//	maxcvertex=currentvertex;      //*
		//}

		for(i=1;i<=gvertexnum;i++)
			if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Coverage[i]>maxcoverage)  //�ֲ�̰�ĵ�ѡ��һ���ڵ�(��֧���)
			{   	maxcoverage=Coverage[i];
					maxcvertex=i;
			}
		if(maxcoverage==0)   //��ʾ������κ��ڵ���Ϊ֧��㶼���������κ�δ����ĵ������̶�. ����֤����������ڲ����������ͼ�в��ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}

		DSList[maxcvertex]=1;          //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //������ص��δ����̶�Satisfied[.]
			if(GAdjMatrix[maxcvertex][j])
			{   Satisfied[j]++;
				if(Satisfied[j]==0)    //�����j��֧����������ܴﵽ����һ��,���������ĵ���
				{	uncoverednum--;
				    for(k=1;k<=gvertexnum;k++)               //������ص�ĸ�����Coverage[.]**************
						if(GAdjMatrix[j][k] && DSList[k]==0)
							Coverage[k]--;
				}
			}
		}


	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewLocalGreedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewLocalGreedy1_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);
	//�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewLocaGreedy1_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewLocalGreedy1_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}


//***����Ƶ��µľֲ�̰�ķ���2***
int NewLocalGreedy2_PIDS(void)
{//����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ���ʹ�ò���2(ʹ��1��).
 //������u��Satisfied[]ֵΪ���Ҿ���ֵ���ĵ�.�ڵ�u������ѡ����Ҫ�̶����(��Needed[.]����ֵ�������u����)��һ����֧���w��Ϊ֧���.
 //����:u����Ҫ������wҲ�����Ҫ
	int i,j,k,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                       //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
	}
    for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		Needed[i]=0;         //***��¼ÿ����֧���(��DSList[.]ֵΪ0�ĵ�)������δ�����ڵ��δ�����֮��(Ϊ����,�����ֵ��ʾ�÷�֧��㱻��Ҫ�ĳ̶�)
		for(j=1;j<=gvertexnum;j++)
           if(GAdjMatrix[i][j])
			   Needed[i]+=Satisfied[j];

	}
    dssize=0;//֧��������ʼ��

	//�ֲ�̰���㷨
	int uncoverednum;
	int currentvertex, minsatisfied, minneeded, minneededvertex;

	uncoverednum=gvertexnum;         //��ʼ��,��ʾ��ʼʱδ�����ĸ���
	while(uncoverednum)
	{//���uncoverednum>0, ������δ����Ľڵ����Ҫ��������µ�֧���

        minsatisfied=0;              //�����ҳ�����̶���С�Ľڵ�(��Satisfied[]ֵΪ����ֵ��С�ĵ�,�����ֵ��ʾ��δ����̶ȵ����)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //Ѱ�ҵ�ǰ����̶���С�Ľڵ㴦��(����򵥵㴦��,Ҳ���԰��������δ����ڵ����ȴ���)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******�ֲ�̰�ķ�����:����ǰδ����Ľڵ�currentvertex�����б���Ҫ�̶�����һ����֧����Ϊ֧���
		minneeded=0;   //�����ҳ�����Ҫ�̶����(��Needed[]ֵΪ���Ҿ���ֵ���)��һ����֧���(currentvertex���ڵ�)

		//����ȡ�����漸�б��Ϊ  //**�Ĵ���(���ʵ�鷢�ֲ����ǵ�ǰ��Ҳ�ɱ�Ϊ֧���ʱЧ�����õĻ�)
		//if(DSList[currentvertex]==0 && Needed[currentvertex]<minneeded)  //** �����ǰδ����ڵ�currentvertex�Ƿ�֧���ʱ,Ҳ�������ɱ�Ϊ֧���
		//{	minneeded=Needed[currentvertex];  //*
		//	minneededvertex=currentvertex;      //*
		//}
		for(i=1;i<=gvertexnum;i++)
			if(GAdjMatrix[currentvertex][i] && DSList[i]==0 && Needed[i]<minneeded)  //�ֲ�̰�ĵ�ѡ��һ���ڵ�(��֧���)
			{   	minneeded=Needed[i];
					minneededvertex=i;
			}
	    if(minneeded==0)   //��ʾ������κ��ڵ���Ϊ֧��㶼���������κ�δ����ĵ������̶�. ����֤����������ڲ����������ͼ�в��ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}

		DSList[minneededvertex]=1;     //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //������ص��δ����̶�Satisfied[.]
			if(GAdjMatrix[minneededvertex][j])
			{    Satisfied[j]++;
 				 if(Satisfied[j]<=0)   //ֻ�ж�֮ǰ��δ������ڵ�(�ոջ������ĵ�֮ǰҲ��δ�����)���б�Ҫ���������Ϣ
				 {	for(k=1;k<=gvertexnum;k++)                //������ص���ڵ�Ĳ�����֮��Needed[]
    					if(GAdjMatrix[j][k] && DSList[k]==0)  //�ڵ��Ƿ�֧���ʱ���б�Ҫ����
							Needed[k]++;                      //ע: ֵΪ��,�����ֵ��ʾ����Ҫ�ĳ̶�
				}
				if(Satisfied[j]==0)                           //�����j��֧����������ܴﵽ����һ��,�����������ĵ���
					 uncoverednum--;
			}
	}

    //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewLocalGreedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewLocalGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);

    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewLocaGreedy2_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewLocalGreedy2_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}


//����Ƶ��µľֲ�̰�ķ���3(�����㷨--���㷨ʱ�临�Ӷ�ΪO(n^2),�����㷨ʱ�临�Ӷ���O(n^3))
int NewLocalGreedy3_PIDS(void)
{//����: ÿ��ѡ��ǰ������u, Ȼ����u�ķ�֧������ڵ��������ȵĵ�Ϊ��֧���.(��ͼ��ͼʱԤ������ͼ�����Ѱ���Ƚ�������)
 //������u��Satisfied[.]ֵΪ���Ҿ���ֵ���ĵ�.��u��������Ӷ�����һ����֧�����ڵ���Ϊ�µ�֧���
	int i,j,dssize;

	for(i=1;i<=gvertexnum;i++)
	{//��ʼ��
		DSList[i]=0;                       //***��¼֧���(ֵΪ0�ĵ�Ϊ��֧���,ֵΪ1�ĵ�Ϊ֧���)
		Satisfied[i]=-(DegreeList[i]+1)/2; //***��¼ÿ�����δ����̶�(ֵ>=0��ʾ���������; ֵ<0��ʾδ�����,�������ֵ��ʾ��δ����ĳ̶�)
	}

    dssize=0;//֧��������ʼ��

	//�ֲ�̰���㷨
	int uncoverednum;
	int currentvertex, minsatisfied, bestvertex;

	uncoverednum=gvertexnum;         //��ʼ��,��ʾ��ʼʱδ�����ĸ���
	while(uncoverednum)
	{//���uncoverednum>0, ������δ����Ľڵ����Ҫ��������µ�֧���

        minsatisfied=0;              //�����ҳ�����̶���С�Ľڵ�(��Satisfied[]ֵΪ����ֵ��С�ĵ�,�����ֵ��ʾ��δ����̶ȵ����)
		currentvertex=0;
		for(i=1;i<=gvertexnum;i++)   //Ѱ�ҵ�ǰ����̶���С�Ľڵ㴦��(����򵥵㴦��,Ҳ���԰��������δ����ڵ����ȴ���)
		    if(Satisfied[i]<minsatisfied)
			{
		    	minsatisfied=Satisfied[i];
		    	currentvertex=i;
			}

		if(currentvertex==0) break;

        //******�ֲ�̰�Ĳ���:����ǰδ����Ľڵ�currentvertex�����е������һ����֧����Ϊ֧���
		for(i=1;i<=gvertexnum;i++)
			if(GAdjMatrix[currentvertex][i] && DSList[i]==0)  //�ֲ�̰�ĵ�ѡ��һ���ڵ�(��֧���)--��ȴ�����(��Ȼ���򼴱�ʾ�ȴ������,������Ϊͼ�����Ѱ���Ƚ�������)
			{		bestvertex=i;
			        break;
			}
		if(i==gvertexnum+1)   //��ʾ������κ��ڵ���Ϊ֧��㶼���������κ�δ����ĵ������̶�. ����֤����������ڲ����������ͼ�в��ᷢ��
		{
			printf("No solution to the instance!\n");
			return 0;
		}

		DSList[bestvertex]=1;          //�����µ�֧���
		dssize++;                      //֧����������
		for(j=1;j<=gvertexnum;j++)     //������ص��δ����̶�Satisfied[.]
			if(GAdjMatrix[bestvertex][j])
			{    Satisfied[j]++;
  				 if(Satisfied[j]==0)   //�����j��֧����������ܴﵽ����һ��,�����������ĵ���
					 uncoverednum--;
			}
	}


    //***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(dssize))
	{	printf("There is some wrong with NewLocalGreedy3_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewLocalGreedy3_PIDS������Ӱ��֧�伯��СΪ(�Ż�ǰ): %5d\n",dssize);

    //�Ż��׶�: ��������Ӱ��֧�伯����һ����С��Ӱ��֧�伯
    int newdssize;
	newdssize=RefinePIDS(dssize);

	//***���ڵ��Խ׶�: ���������Ӱ���֧�伯�Ƿ���ȷ
	if(!CheckPIDS(newdssize))
	{	printf("There is some wrong with NewLocaGreedy3_PIDS!\n");
		getch();
	}//***���ڵ��Խ׶�

	printf("\nNewLocalGreedy3_PIDS������Ӱ��֧�伯��СΪ(�Ż���): %5d\n",newdssize);
	return newdssize;//��������֧�伯��С
}







//***������������***

void ReadGraph(char *txtgraphfilename)
{//To read the data file of graph G=<V,E>(V={1,2,...,n}) to obtain the adjacent matrix.
 //The format of the data file is as follows: the node number,the edge number, every edge <x,y> followed
      FILE  *in;

      //���ݸ������ļ����ִ��ļ�������r��ʾ���ļ�����ֻ������
	  if ((in=fopen(txtgraphfilename, "r"))==NULL)
      {  	fprintf(stderr, "Cannot open the data file.\n");
		return;
      }

      fscanf(in,"%d",&gvertexnum);                         //��ȡͼ��ʵ�ʵ���
      fscanf(in,"%d",&gedgenum);                           //��ȡͼ�ı���

      int i,j,k;
      for(i=1;i<=gvertexnum;i++)                          //initialize the adjacent matrix
		for(j=1;j<=gvertexnum;j++)
			GAdjMatrix[i][j]=0;             /**/

      k=0;                                                 //k ������סʵ�ʵı�����ע���ļ���ĳ�߿����ظ����ֶ��
      while(!feof(in))                                     //����ÿһ����(ÿ�����ļ�������һ�Ե�����ʾ��)
      {
	    fscanf(in,"%d",&i);
    	fscanf(in,"%d",&j);

	    if(GAdjMatrix[i][j])                              //��������˾Ͳ�Ҫ�ظ���ע���ļ��б߿����ظ����ֶ��
             continue;
	    GAdjMatrix[i][j]=1;
	    GAdjMatrix[j][i]=GAdjMatrix[i][j];                //For an undirected graph, the adjacent matrix is symmetrical
	    k++;
      }
      gedgenum=k;                                         //�õ�ʵ�ʵı���
      fclose(in);                                         //�򿪵��ļ�������ر�

     //��ʼ�������,�����������(ͼ�к��й������������޽�). ����Ƿ��й�����. ������й�����,�򽫹�����������߼��ɱ�֤����������
  	 for(i=1;i<=gvertexnum;i++)
	 {  DegreeList[i]=0;
		for(j=1;j<=gvertexnum;j++)
			if(GAdjMatrix[i][j])
				DegreeList[i]++;
		if(DegreeList[i]==0)  //��ʾ�ڵ�i�ǹ�����,��ʱ���ڵ�i������ӵ���һ�ڵ�k
		{
           k=(rand()%gvertexnum)+1;
		   if(k==i && k<gvertexnum)
			   k++;
		   else if(k==i && k>1)
			   k--;
           GAdjMatrix[i][k]=1;
		   GAdjMatrix[k][i]=1;
		   DegreeList[i]++;
		   DegreeList[k]++;
		   gedgenum++;
		}
	 }

	 //�Խڵ�Ȱ����ɴ�С����, Ȼ���ͼ���ն��ɴ�С���±��: ��1�Žڵ�����, ������gvertexnum�Žڵ����С
	 //����ѡ�����򷽷�
	 short int temp;
	 bool tempf;
     int a;
	 for(i=1;i<gvertexnum;i++)
	 {
		 k=i;
		 for(j=i+1;j<=gvertexnum;j++)
			 if(DegreeList[j]>DegreeList[k])
			     k=j;
	     if(k!=i)
		 {//�����ڵ�i�ͽڵ�k�ı��: ��Ӧ�ڵ��Ҫ����,ͬʱҪ�����ڽӾ���ĵ�i�к͵�k��,�Լ���i�к͵�k��
            temp=DegreeList[k];        //�������ڵ�Ķ�
			DegreeList[k]=DegreeList[i];
			DegreeList[i]=temp;

			for(a=1;a<=gvertexnum;a++) //�����ڽӾ��������
			{	tempf=GAdjMatrix[k][a];
				GAdjMatrix[k][a]=GAdjMatrix[i][a];
				GAdjMatrix[i][a]=tempf;
			}

            for(a=1;a<=gvertexnum;a++) //�����ڽӾ��������
			{	tempf=GAdjMatrix[a][k];
				GAdjMatrix[a][k]=GAdjMatrix[a][i];
				GAdjMatrix[a][i]=tempf;
			}

		 }
	 }

	 //������maxdegree, ��С��mindegree, ���ж�Ϊ1�Ľڵ����onedegreenum, ż���Ƚڵ����evendegreenum(�����Ƚڵ����Ϊgvertexnum-evendegreenum)
	 maxdegree=DegreeList[1];
	 mindegree=DegreeList[gvertexnum];
	 onedegreenum=0;
	 evendegreenum=0;
	 for(i=gvertexnum;i>=1;i--)
		 if(DegreeList[i]==1)
			 onedegreenum++;
		 else if(DegreeList[i]%2==0)
			 evendegreenum++;

	 printf("n=%d,m=%d,maxd=%d,averd=%6.2lf,dens=%10.5lf\n",gvertexnum,gedgenum,maxdegree,2.0*gedgenum/gvertexnum,2.0*gedgenum/(gvertexnum*(gvertexnum-1.0)));
	 printf("mindegree=%d,onedegreenum=%d,evendegreenum=%d,odddegreenum=%d\n",mindegree,onedegreenum,evendegreenum,gvertexnum-evendegreenum);

}





//���ɻ�����ִ�в���(begin)****************************************************************
  int main(void)
  {

    char  graphfilename[50];
	char  resultfilename[50];
    //������������ReadGraph()��Ҫ�򿪵�ͼ�����ļ���(txt�ļ�,����չ��.�ļ���ʽ:��һ��ͼ�ĵ���n,�ڶ���ͼ�ı���m,����ÿһ����һ����(i,j)��һ���ڵ�i  j)
	printf("Please input the filename of the data file:\n");
	scanf("%s",graphfilename);
    //������������������Ҫ�����ı�����������ļ���(txt�ļ�,����չ��)
	printf("Please input the filename of the result file:\n");
	scanf("%s",resultfilename);

//���ɻ�����ִ�в���(end)******************************************************************





	//------------------------------��������--------------------------------------

  	ReadGraph(graphfilename);  //����ͼ������
    //CreateGraph(graphfilename);  //�������һ��ͼ,���洢���ļ���


	//-----------------------------ִ���㷨--------------------------------------
	FILE *out;
	clock_t start, end;
    int dssize;

	//���ݸ����Ľ���ļ����ִ��ļ�������w��ʾ���ļ�����д����
    if ((out=fopen(resultfilename, "w"))==NULL)
	{   	fprintf(stderr, "Cannot open the data file.\n");
			return 0;
	}
	fprintf(out,"******The graph file is %s******\n",graphfilename);
	fprintf(out,"vertexnum=%d,edgenum=%d,maxdegree=%d,averdegree=%6.2lf,density=%10.5lf\n",gvertexnum,gedgenum,maxdegree,2.0*gedgenum/gvertexnum,2.0*gedgenum/(gvertexnum*(gvertexnum-1.0)));
	fprintf(out, "mindegree=%d,onedegreenum=%d,evendegreenum=%d,odddegreenum=%d\n\n",mindegree,onedegreenum,evendegreenum,gvertexnum-evendegreenum);

    //******̰�Ľ����㷨1(�����з���)*********
	printf("\nGreedy1_PIDS is running...");
    start=clock();
   	dssize=Greedy1_PIDS();      //ִ��̰�Ľ����㷨,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
	fprintf(out,"\nGreedy1_PIDS has running time %.3lf.\n", 1.0*(end-start)/CLK_TCK);
	//getch();

	/*
	//******̰���㷨2(�����з���)*********
	printf("\nGreedy2_PIDS is running...");
    start=clock();
   	dssize=Greedy2_PIDS();      //ִ��̰�Ľ����㷨,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
	fprintf(out,"\nGreedy2_PIDS has running time %.3lf.\n", 1.0*(end-start)/CLK_TCK);
	//getch();


	//*****����Ƶ�̰���㷨1*********
	printf("\nNewGreedy1_PIDS is running...");
    start=clock();
	dssize=NewGreedy1_PIDS();      //ִ������Ƶ�̰���㷨1,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nNewGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewGreedy1_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶ�̰���㷨2*********
	printf("\nNewGreedy2_PIDS is running...");
    start=clock();
	dssize=NewGreedy2_PIDS();      //ִ������Ƶ�̰���㷨2,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nNewGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewGreedy2_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶ�̰���㷨2*********
	printf("\nNewGreedy3_PIDS is running...");
    start=clock();
	dssize=NewGreedy3_PIDS();      //ִ������Ƶ�̰���㷨2,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nNewGreedy3_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewGreedy3_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶľֲ�̰���㷨1*********
	printf("\nLocalGreedy1_PIDS is running...");
    start=clock();
	dssize=LocalGreedy1_PIDS();      //ִ������Ƶľֲ�̰���㷨1,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nLocalGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nLocalGreedy1_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶľֲ�̰���㷨2*********
	printf("\nLocalGreedy2_PIDS is running...");
    start=clock();
	dssize=LocalGreedy2_PIDS();      //ִ������Ƶľֲ�̰���㷨2,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nLocalGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nLocalGreedy2_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶľֲ�̰���㷨3*********
	printf("\nLocalGreedy3_PIDS is running...");
    start=clock();
	dssize=LocalGreedy3_PIDS();      //ִ������Ƶľֲ�̰���㷨3,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nLocalGreedy3_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nLocalGreedy3_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶ��µľֲ�̰���㷨1*********
	printf("\nNewLocalGreedy1_PIDS is running...");
    start=clock();
	dssize=NewLocalGreedy1_PIDS();  //ִ������Ƶ��µľֲ�̰���㷨1,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nNewLocalGreedy1_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewLocalGreedy1_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶ��µľֲ�̰���㷨2*********
	printf("\nNewLocalGreedy2_PIDS is running...");
    start=clock();
	dssize=NewLocalGreedy2_PIDS();  //ִ������Ƶ��µľֲ�̰���㷨2,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nNewLocalGreedy2_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewLocalGreedy2_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();

	//*****����Ƶ��µľֲ�̰���㷨3*********
	printf("\nNewLocalGreedy3_PIDS is running...");
    start=clock();
	dssize=NewLocalGreedy3_PIDS();  //ִ������Ƶ��µľֲ�̰���㷨3,��������֧�伯��С
	end=clock();
	if(dssize)
		fprintf(out,"\nNewLocalGreedy3_PIDS obtains a PIDS of size %d.", dssize);
	else
		fprintf(out,"\nNo solution to the instance.");
    fprintf(out,"\nNewLocalGreedy3_PIDS has running time %.3lf\n", 1.0*(end-start)/CLK_TCK);
    //getch();


    */

	fclose(out);           //�򿪵��ļ�������ر�
    printf("�����������!\n");
	//getch();
}































/*�������ͼ�ĺ���û��ʹ�á���������*/
void CreateGraph(char *txtgraphfilename)
{//��������ͼ,ͼ��n����1,2,...,n, ͼ��ƽ����averd=2m/n, �ߵĳ��ܶ�r=2m/(n(n-1)),�����r=averd/(n-1).
  //ͼ�������i��j֮����ڱߵĸ���Ϊr(0<r<1)--r���Կ����Ǳߵĳ��ܶ�.
   printf("Please input the number of nodes of the graph:\n");
   scanf("%d",&gvertexnum);       //ȷ��ͼ�ĵ���

   double r;
   int averd;
   //printf("Please input the density r of edges of the graph (0<r<1):\n");
   //scanf("%lf",&r);             //ȷ��ͼ�бߵĳ��ܶ�
   printf("Please input the average degree of the graph (1<=averd<=n-1):\n");
	   scanf("%d", &averd);       //����ͼ��ƽ�����
   r=(1.0*averd)/(gvertexnum-1.0);           //ȷ��ͼ�бߵĳ��ܶ�

   int i,j;
   for(i=1;i<=gvertexnum;i++)     //��ʼ���ڽӾ���
	   for(j=1;j<=gvertexnum;j++)
  	      GAdjMatrix[i][j]=0;

   time_t t;

   //�������ͼ. ������й������򽫹�����������߼����Ա�֤����������

   srand((unsigned) time(&t));
   gedgenum=0;

   for(i=1;i<gvertexnum;i++)     //�������ͼ���ڽӾ���
	   for(j=i+1;j<=gvertexnum;j++)
	   {
		   GAdjMatrix[i][j]=rand()%1000<(r*1000)?1:0;
	       if(GAdjMatrix[i][j])
		   {
			   GAdjMatrix[j][i]=1;
		       gedgenum++;
		   }
	   }

//��ʼ�������,�����������(ͼ�к��й������������޽�). ����Ƿ��й�����. ������й�����,�򽫹�����������߼��ɱ�֤����������
	 int k;
  	 for(i=1;i<=gvertexnum;i++)
	 {  DegreeList[i]=0;
		for(j=1;j<=gvertexnum;j++)
			if(GAdjMatrix[i][j])
				DegreeList[i]++;
		if(DegreeList[i]==0)  //��ʾ�ڵ�i�ǹ�����,��ʱ���ڵ�i������ӵ���һ�ڵ�k
		{
           k=(rand()%gvertexnum)+1;
		   if(k==i && k<gvertexnum)
			   k++;
		   else if(k==i && k>1)
			   k--;
           GAdjMatrix[i][k]=1;
		   GAdjMatrix[k][i]=1;
		   DegreeList[i]++;
		   DegreeList[k]++;
		   gedgenum++;
		}
	 }

	 //�Խڵ�Ȱ����ɴ�С����, Ȼ���ͼ���ն��ɴ�С���±��: ��1�Žڵ�����, ������gvertexnum�Žڵ����С
	 //����ѡ�����򷽷�
	 short int temp;
	 bool tempf;
     int a;
	 for(i=1;i<gvertexnum;i++)
	 {
		 k=i;
		 for(j=i+1;j<=gvertexnum;j++)
			 if(DegreeList[j]>DegreeList[k])
			     k=j;
	     if(k!=i)
		 {//�����ڵ�i�ͽڵ�k�ı��: ��Ӧ�ڵ��Ҫ����,ͬʱҪ�����ڽӾ���ĵ�i�к͵�k��,�Լ���i�к͵�k��
            temp=DegreeList[k];        //�������ڵ�Ķ�
			DegreeList[k]=DegreeList[i];
			DegreeList[i]=temp;

			for(a=1;a<=gvertexnum;a++) //�����ڽӾ��������
			{	tempf=GAdjMatrix[k][a];
				GAdjMatrix[k][a]=GAdjMatrix[i][a];
				GAdjMatrix[i][a]=tempf;
			}

            for(a=1;a<=gvertexnum;a++) //�����ڽӾ��������
			{	tempf=GAdjMatrix[a][k];
				GAdjMatrix[a][k]=GAdjMatrix[a][i];
				GAdjMatrix[a][i]=tempf;
			}

		 }
	 }

	 //������maxdegree, ��С��mindegree, ���ж�Ϊ1�Ľڵ����onedegreenum, ż���Ƚڵ����evendegreenum(�����Ƚڵ����Ϊgvertexnum-evendegreenum)
	 maxdegree=DegreeList[1];
	 mindegree=DegreeList[gvertexnum];
	 onedegreenum=0;
	 evendegreenum=0;
	 for(i=gvertexnum;i>=1;i--)
		 if(DegreeList[i]==1)
			 onedegreenum++;
		 else if(DegreeList[i]%2==0)
			 evendegreenum++;

	 printf("n=%d,m=%d,maxd=%d,averd=%6.2lf,dens=%10.5lf\n",gvertexnum,gedgenum,maxdegree,2.0*gedgenum/gvertexnum,2.0*gedgenum/(gvertexnum*(gvertexnum-1.0)));
	 printf("mindegree=%d,onedegreenum=%d,evendegreenum=%d,odddegreenum=%d\n",mindegree,onedegreenum,evendegreenum,gvertexnum-evendegreenum);

	bool savegraph;
	savegraph=1;                      //savegraph=0 ��ʾ���ô洢ͼ���ļ���
	if(savegraph)
	{	//�洢ͼ��ָ�����ļ���
		FILE  *out;
		//���ݺ���������ļ����֣����������ļ�������w��ʾ���ļ�����д����
		if ((out=fopen(txtgraphfilename, "w"))==NULL)
		{   	fprintf(stderr, "Cannot open the data file.\n");
			return;
		}
		fprintf(out,"%d\n",gvertexnum);
		fprintf(out,"%d\n",gedgenum);

		for(i=1;i<gvertexnum;i++)
			for(j=i+1;j<=gvertexnum;j++)
				if(GAdjMatrix[i][j])    //�洢ÿ����
					fprintf(out,"%5d  %5d\n",i,j);
		fclose(out);                   //�򿪵��ļ�������ر�
	}
}


