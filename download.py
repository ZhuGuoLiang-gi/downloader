import os
import requests
from tqdm import tqdm
import time
import threading
import sys
import multiprocessing
import argparse
import psutil
import time
import datetime
from rich.progress import Progress
import argparse



'''
this script description and information:
author:zhuguoliang
Email:zhuguoliang001@gmail.com
'''



'''
bioinformation database about protein、DNA、RNA:
structure or sequences information website:
uniprot website:https://www.uniprot.org/help/downloads
    Embeddings:
        Reviewed (Swiss-Prot)	https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/uniprot_sprot/per-protein.h5
        Arabidopsis thaliana	https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000006548_3702/per-protein.h5
        Caenorhabditis elegans	https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000001940_6239/per-protein.h5
        Escherichia coli	    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000000625_83333/per-protein.h5
        Homo sapiens	        https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000005640_9606/per-protein.h5
        Mus musculus	        https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000000589_10090/per-protein.h5
        Rattus norvegicus	    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000002494_10116/per-protein.h5
        SARS-CoV-2	            https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000464024_2697049/per-protein.h5
        
    UniRef
        UniRef100               https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
        UniRef90                https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
        UniRef50                https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
'''



class File():
    def __init__(self,name,size,url='',path='./'):
        self.name=name
        self.size=size
        self.star=0
        self.end=size-1
        self.child=[]
        self.url=url
        self.path=path
        self.file=os.path.join(path,name)
        self.parent=''
    
    def p_n(self,memory,pro_num):
        pro_num_max = self.size // memory
        if pro_num_max == 0:
            return 1
        if pro_num_max > pro_num:
            return pro_num
        if pro_num_max < pro_num:
            return pro_num_max
        return pro_num
        
    def th_n(self,thread_number):
        min_menory = 1024 *1024 * 3
        max_thread = 50
        if min_menory >= self.size // thread_number:
            return 1
        if min_menory < self.size // thread_number and thread_number < max_thread :
            return thread_number
        if thread_number > max_thread :
            return max_thread
        return thread_number
        
    def split_N(self,N):
        if N == 1:return []
        child_file=[]
        space_size=self.size // N
        for i in range(N):
            child_file.append(File('',0))
            child_file[-1].name=f'child_{i}_{self.name}'
            child_file[-1].star=i*space_size + self.star
            child_file[-1].end=(i+1)*space_size-1 + self.star if i < N -1 else self.end
            child_file[-1].size=child_file[-1].end-child_file[-1].star + 1
            child_file[-1].url=self.url
            child_file[-1].path=self.path
            child_file[-1].parent=self
            child_file[-1].file=os.path.join(self.path,child_file[-1].name)
            self.child.append(child_file[-1])
        return child_file
    
    def check(self): 
        
        # 1:file is ok
        # 0:file has problem
        
        if os.path.exists(self.file):
            if os.path.getsize(self.file) >= self.size:
                return 1
            else:
                return 0
        else:
            if not self.parent:return 0
            if os.path.exists(self.parent.file):
                if os.path.getsize(self.parent.file) == self.parent.size:
                    return 1
                else:
                    return 0
            else:
                return 0
    
    def check_child(self):
        
        # 0:child file of downloading has problem 
        # 1:child file of downloading file is ok
        # 2:downloading file is node  (not exist child file)
        
        if len(self.child)==0:
            return 2
        for child in self.child:
            if child.check()==0:
                return 0
        return 1
    
    def merge_child(self):
        with open(self.file,'wb') as parent_file:
            for child in self.child:
                with open(child.file,'rb') as child_file:
                    for line in child_file.readlines():
                        parent_file.write(line)
        print()
        print(f'{self.name} merge successfully')
        print(f"[info] {datetime.datetime.now().replace(microsecond=0)}: delete child file of {self.name}")
        print()
        self.delect_child()
        return
    
    def delect_child(self):
        for child in self.child:
            os.remove(child.file)
        return
    
    def download(self,progress_bar=None):
        
        if not progress_bar:
            progress = Progress()
        else:
            progress=progress_bar
        
        headers = {'Range': f'bytes={self.star}-{self.end}'}
        retry_delay = 3
        try_max=5
        try_num=3
        

        while True:
            try:
                response = requests.get(self.url,headers=headers, stream=True)
                break
            except:
                time.sleep(retry_delay)
                try_num+=1
                if try_num > try_max:
                    print(f'{self.name} download failed !!!')
                    exit()
                return
            
        chunk_size = 1024 * 2
          
        with progress:
            
            task = progress.add_task(f"[cyan] {self.name} download ...", total=self.size)
            completed=0
        
            chunk_try=0
            chunk_try_max=10
            with open(self.file, 'wb') as f:
                
                response_star=time.time()
                response_time=0.02
                response_max_time=3
                response_times=0
            
                for chunk in response.iter_content(chunk_size=chunk_size):
                
                    response_end=time.time()
                
                    if response_end - response_star > response_time:
                        response_limit_time+=1
                        if response_times > response_max_time:
                            print(f'[warnning] {self.name} response time is too long !!! ')
                            return 0
                    else:
                        response_limit_time=0
                
                    if chunk:
                        f.write(chunk)
                        completed=completed+len(chunk)
                        progress.update(task, completed=completed,description=f"{self.name} is downloading... ({completed / 1024 / 1024:.2f}/{self.size / 1024 / 1024:.2f} MB)")
                    else:
                        chunk_try+=1
                        print(f'{self.name},test 2')
                        if chunk_try > chunk_try_max:
                            print("download has problem !!! ....................................")
                            return 0
                    response_star=time.time()
        return
    
    def run_download(self,thread_num=1):    
                  
        if thread_num==1:
            self.download()
            return         
        
        if thread_num > 1:
            
            progress=Progress()
            
            self.split_N(thread_num)
            threads=[]
            fail=[]
            for child in self.child:
                if not child.check(): fail.append(child)
                    
            if len(fail) == 0:
                print(f'{self.name} download successfully')
            elif len(fail) == 1:
                fail[-1].download()
            elif len(fail)>1:
                for child in fail:
                    thread=threading.Thread(target=File.download,args=(child,progress))
                    threads.append(thread)
                
                for thread in threads:
                    thread.start()
                    
                for thread in threads:
                    thread.join()
        return
        
    def run_merge(self):
        if self.check() == 1:
            print(f'{self.name} has merge successfully')
        else:
            if self.check_child()==1:
                self.merge_child()
                print(f'{self.name} merge successfully')
            if self.check_child()==2:
                print(f'{self.name} merge fail,need to download again')
            if self.check_child()==0:
                print(f'child file of {self.name} has problem')
        return
                
    def solve(self):
        if self.check():
            return 1
        else:
            if self.check_child()==1:pass
            if self.check_child()==2:pass
            if self.check_child()==0:
                if os.path.exists(self.file):os.remove(self.file)
                for child in self.child:
                    if os.path.exists(child.file):
                        if child.check_child() == 1:pass
                        if child.check_child() == 2:pass
                        if child.check_child() == 0:
                            if os.path.exists(child.file):os.remove(child.file)
                            for grand_child in child.child:
                                pass
                                #if os.path.exists(os.path.join(grand_child.path,grand_child.name)): os.remove(os.path.join(grand_child.path,grand_child.name))
        return 0
    

def download_file(child,thrumd_num):
    child.run_download(thrumd_num)
    child.run_merge()
    return 0

def main(file,pro_num,thread_num):
    if file.child:
        with multiprocessing.Pool(pro_num) as pool:
            task=[(child,thread_num) for child in file.child]
            Pool = pool.starmap(download_file,task)
    # Pool.get()
    if not file.child:
        download_file(file,thread_num)
        
    file.run_merge()
    return 0

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="this is a download tool to download huge file faster")
    
    parser.add_argument('-u','--url',dest='url',required=True, help='input downloaded file address')
    parser.add_argument('-o','--out',dest='out',default="./",help='download file to path')
    parser.add_argument('-p','--pro_num',dest='pro_num',default=4,help='set processes number')
    parser.add_argument('-t','--thr_num',dest='thr_num',default=10,help='set thread number')

    args = parser.parse_args()

    url = args.url
    pro_num=int(args.pro_num)
    thread_num=int(args.thr_num)
    path=args.out
    menory=1024*1024*100

    start_time=time.time()
    response = requests.head(url)
    fild_name = url.split("/")[-1]
    try:
        file_size = int(response.headers['Content-Length'])
    except:
        file_size= 1024
 
    
    file=File(fild_name,file_size,url,path)
    pro_num=file.p_n(menory,pro_num)
    thread_num=file.th_n(thread_num)
    file.split_N(pro_num)
    
    try_max=5
    try_num=0
    
    while not file.solve():
        try_num+=1
        print(f'try to download:{try_num}')
        
        if try_num > try_max:
            print(f'=================================================')
            print('')
            print(f'downloading file have some problem!!')
            print(f'you could executate this program aganin to fix it')
            print('')
            print(f'==================================================')
            exit()
        
        main(file,pro_num,thread_num)
    
    end_time = time.time()
    execution_time = end_time - start_time
    min=execution_time / 60
    print(f'running time: {min:.2f} min')
    



        


