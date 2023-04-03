import requests
import jsonpath
import argparse
import os


class Gcat:
    def __init__(self, qid_path="", out_dir=""):
        self.qid = ""
        self.qidpath = qid_path
        self.out_dir = out_dir
        self.tName = ""
        self.urllist = []
        self.outlist = []
        self.headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"}
        
    def get_traitN(self):
        data = {"q": self.qid}
        res = requests.post("https://www.ebi.ac.uk/gwas/api/search/advancefilter", headers=self.headers, data=data)
        page_dict = res.json()
        self.tName = jsonpath.jsonpath(page_dict, "$..traitName")[0][0]
        
        
    def download(self):
        url_startnum = int(self.qid.replace("GCST", ""))//1000 * 1000 +1
        url_endnum = url_startnum + 999
        url1 = f"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{url_startnum}-GCST{url_endnum}/{self.qid}/{self.qid}_buildGRCh38.tsv.gz"
        url2 = f"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{url_startnum}-GCST{url_endnum}/{self.qid}/{self.qid}_buildGRCh38.tsv.gz"
        url3 = f"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{url_startnum}-GCST{url_endnum}/{self.qid}/{self.qid}_buildGRCh37.tsv"
        url4 = f"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{url_startnum}-GCST{url_endnum}/{self.qid}/{self.qid}_buildGRCh38.tsv"
        self.urllist = [url1, url2, url3, url4]
        output1_basename = (self.qid+"_"+self.tName+"_buildGRCh38.tsv.gz").replace(" ", "_").replace("(", "_").replace(")", "_")
        output1 = os.path.join(self.out_dir, output1_basename)
        output2_basename = (self.qid+"_"+self.tName+"_buildGRCh37.tsv.gz").replace(" ", "_").replace("(", "_").replace(")", "_")
        output2 = os.path.join(self.out_dir, output2_basename)
        output3_basename = (self.qid+"_"+self.tName+"_buildGRCh38.tsv").replace(" ", "_").replace("(", "_").replace(")", "_")
        output3 = os.path.join(self.out_dir, output3_basename)
        output4_basename = (self.qid+"_"+self.tName+"_buildGRCh37.tsv").replace(" ", "_").replace("(", "_").replace(")", "_")
        output4 = os.path.join(self.out_dir, output4_basename)
        self.outlist = [output1, output2, output3, output4]
        flag = False
        for output, url in zip(self.outlist, self.urllist):
            scode = requests.head(url, headers=self.headers).status_code
            if scode == 200:
                flag = True
                print(f"find valid {url} for trait {self.tName}")
                os.system(f"wget -O {output} -c {url}")
        if not flag:
            print(f"Not found url for {self.qid} of {self.tName} ")
            
    def pipeline(self):
        self.get_traitN()
        self.download()
    
    def loop_pipeline(self):
        if self.qidpath:
            with open(self.qidpath, "r") as fp:
                for id in fp:
                    self.qid = id.rstrip()
                    self.pipeline()
        print("All finished!")

def getcmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--idlist", type=str, help="Path to a file include GCST_XXX ids")
    parser.add_argument("-o", "--output", type=str, help="Specfiy output directory", default="./")
    args = parser.parse_args()
    qid_path = args.idlist
    out_dir = args.output
    return qid_path, out_dir


if __name__ == "__main__":
    qid_path, out_dir = getcmd()
    gcat = Gcat(qid_path, out_dir)
    gcat.loop_pipeline()
