"""Render the standalone CSA dashboard and exercise its real browser controls."""
from pathlib import Path
import argparse,html,json,re,subprocess
ROOT=Path(__file__).resolve().parents[2]

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('html',type=Path);a=p.parse_args()
 work=ROOT/'notes/physics/work/ion_csa_html_test';work.mkdir(parents=True,exist_ok=True)
 javascript=r'''<script>
window.addEventListener('error',e=>{document.body.dataset.testError=String(e.message)});
(async()=>{let checks=0;const pause=()=>new Promise(r=>setTimeout(r,150));try{
 await pause();
 for(let id of ['snr','eventsnr','stories','bias','cov','noise','story']){if(!document.getElementById(id).data?.length)throw Error('Empty plot '+id);checks++}
 const select=(id,value)=>{let e=document.getElementById(id);e.value=value;e.dispatchEvent(new Event('change'))};
 let previousCards=document.getElementById('cards').textContent;select('sample',groupKeys[1]);await pause();
 if(document.getElementById('cards').textContent===previousCards)throw Error('Sample metrics did not change');checks++;
 select('severity','0');await pause();if(document.getElementById('stories').data[0].type!=='bar')throw Error('Discrete story plot');checks++;
 select('snrscale','low');await pause();if(document.getElementById('snr').layout.xaxis.type==='log')throw Error('Signed S/N view stayed logarithmic');checks++;
 select('normalization','count');await pause();if(document.getElementById('snr').layout.yaxis.title.text!=='Entries per bin')throw Error('Histogram normalization');checks++;
 select('example',exampleKeys[exampleKeys.length-1]);await pause();if(document.getElementById('story').data.length!==5)throw Error('Story traces missing');checks++;
 select('window','all');await pause();let e=DATA.groups[document.getElementById('sample').value].examples[document.getElementById('example').value];if(document.getElementById('story').layout.xaxis.range[1]!==e.edge_time_s[e.edge_time_s.length-1])throw Error('Recovery range');checks++;
 select('sample',groupKeys[0]);select('severity','2');select('snrscale','log');select('normalization','probability');select('window','event');await pause();
 document.body.dataset.testResult=JSON.stringify({passed:true,checks});
}catch(e){document.body.dataset.testResult=JSON.stringify({passed:false,error:String(e),browser_errors:window._errors,checks})}})();
</script>'''
 target=work/'dashboard_test.html';target.write_text(a.html.read_text().replace('<head>',"<head><script>window._errors=[];window.addEventListener('error',e=>window._errors.push(e.message));</script>").replace('</body>',javascript+'</body>'))
 browser=Path('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
 result=subprocess.run([str(browser),'--headless','--disable-gpu','--no-first-run','--allow-file-access-from-files','--user-data-dir='+str(work/'profile'),'--window-size=1440,1100','--virtual-time-budget=6000','--screenshot='+str(work/'dashboard.png'),'--dump-dom',target.as_uri()],capture_output=True,text=True,timeout=90)
 (work/'browser.log').write_text(result.stderr)
 assert result.returncode==0,result.stderr[-1500:]
 match=re.search(r'data-test-result="([^"]+)"',result.stdout)
 assert match,'Dashboard did not complete control tests: '+result.stderr[-1000:]
 status=json.loads(html.unescape(match[1]));assert status['passed'],status
 (work/'result.json').write_text(json.dumps(status,indent=2)+'\n');print(status)
if __name__=='__main__':main()
