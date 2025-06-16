#!/Users/arri/miniconda3/bin/python


# https://python-redmine.com/introduction.html
import redminelib
from redminelib import Redmine

import datetime
import argparse
import re
import textwrap

# import textwrap
# indent = textwrap.TextWrapper(initial_indent='          =  ',subsequent_indent='             ')

# TODOs:
# help
# check for open / closed
# level = 1  ->  date ausgeben
# nur description ausgeben
# selectieren nach priority und status



def get_structure(redmine) :
  # https://python-redmine.com/resources/project.html#all
  #projects = redmine.project.all(limit=1)
  projects = redmine.project.all()
  for p in projects :
    print('-----------------------------------------------------------------------------------------------------------------------------------')
    #print(p.get_parent['id'])
    for i in p :
      print('  %-25s  %-50s  %-30s '%(i[0],type(i[1]),str(i[1])))
      #if i[0]=='parent' :
      #  print(i[1]['name'])
      #if i[1]==None :
      #  for j in i :
      #    print(j)
      #    #print('     %-25s  %-50s  %-30s '%(j[0],type(j[1]),str(j[1])))
  # id                         <class 'int'>                                      
  # name                       <class 'str'>                                      
  # updated_on                 <class 'str'>                                      
  # created_on                 <class 'str'>                                      
  # identifier                 <class 'str'>                                      
  # description                <class 'str'>                                      
  # status                     <class 'int'>                                      
  # is_public                  <class 'bool'>                                     
  # url                        <class 'str'>                                      
  # internal_id                <class 'int'>                                      
  # manager                    <class 'redminelib.managers.standard.ProjectManager'> 
  # wiki_pages                 <class 'NoneType'>                                 
  # memberships                <class 'NoneType'>                                 
  # issue_categories           <class 'NoneType'>                                 
  # time_entries               <class 'NoneType'>                                 
  # versions                   <class 'NoneType'>                                 
  # news                       <class 'NoneType'>                                 
  # issues                     <class 'NoneType'>                                 
  # files                      <class 'NoneType'>                                 
  # trackers                   <class 'NoneType'>                                 
  # enabled_modules            <class 'NoneType'>                                 
  # time_entry_activities      <class 'NoneType'>                                 
  # issue_custom_fields        <class 'NoneType'>                                 



def get_all_projects(redmine) :
  #projects = redmine.project.all(limit=2)
  projects = redmine.project.all()
  projectslist = list(projects)
  #print(projectslist)
  my = 0
  projdata=[]
  for p in projects :
    parent=''
    for i in p :
      if i[0]=='parent' :
        parent=i[1]['name']
    projdata.append( [ p.id, p.name, p.updated_on, parent ] )
  def myKey1(e):
    return e[1]
  projdata.sort(key=myKey1)
  def myKey3(e):
    return e[3]
  projdata.sort(key=myKey3,reverse=True)
  for i in projdata:
    print( '  %-3s  %-40s %16s    %5s ' % (i[0],i[1],i[2],i[3]) )

    

def print_issue(issue,level,last,assigned_to) :

  #print(list(issue.custom_fields.values()))
  #for field in issue.custom_fields.values() :
  #  print(field)

  #print(issue.id)
  #print(issue.done_ratio)

  #if (issue.id==2703) :
  #  redmine.issue.update(issue.id,done_ratio=50.)

  issue_assigned_to = ""
  try:
    issue_assigned_to = str(issue.assigned_to)
  except :
    issue_assigned_to = "--"
    pass

  if assigned_to!='' and assigned_to!=issue_assigned_to :
    return
  
  estimated_hours = ""
  try:
    estimated_hours = issue.estimated_hours
  except :
    estimated_hours = "--"
    pass

  done_ratio = ""
  try:
    done_ratio = issue.done_ratio
  except :
    done_ratio = "--"
    pass
  
  instrument = ''
  objectname = ''
  ra = 0.
  dec = 0.
  airmass = 99.
  lunar_phase = ''
  sky = 99.
  seeing = 99.
  transparency = 99.
  
  if level>=2 :
    print("======================================================================================================================================================")
    print('%-3s  %-7s  %-20s  %-10s  %-11s  %6s  %6s  %7s  %7s  %6s  %14s  %14s  %9s  %4s' % ( 'pID','issueID', 'Objectname', 'Instrument', 'Status', 'RA', 'DEC', 'Sky_min', 'PSF_max', 'AM_max', 'moonphase', 'transparency' ,'est.hours', 'done'))
    print("------------------------------------------------------------------------------------------------------------------------------------------------------")
  for field in issue.custom_fields.values() :

    #print(field)

    if field['name'] == 'Instrument' :
      instrument = field['value'][0]

    if field['name'] == 'Objectname' :
      objectname = field['value']

    if field['name'] == 'RA' :
      ra = field['value']
      if ra != '':
        ra = float(ra)
      else:
        ra = -1.0
        
    if field['name'] == 'DEC' :
      dec = field['value']
      if dec != '':
        dec = float(dec)
      else:
        dec = -1.0

    if field['name'] == 'Airmass (max.)' :
      airmass = field['value']
      if airmass != '':
        airmass = float(airmass)
      else:
        airmass = -1.0

    if field['name'] == 'Lunar Phase' :
      lunar_phase = field['value']

    if field['name'] == 'Sky (min.)' :
      sky = field['value']
      if sky is not None and sky != '' :
        sky = float(sky)
      else:
        sky = -1.0

    if field['name'] == 'Seeing (max.)' :
      seeing = field['value']
      if seeing != '':
        seeing = float(seeing)
      else:
        seeing = -1.0

    if field['name'] == 'Transparency' :
      transparency = field['value']


  if level>=1 :
    print('%-3.0f  %-7.0f  %-20s  %-10s  %-11s  %6.1f  %6.1f  %7.1f  %7.1f  %6.1f  %14s  %14s  %9s  %4s' %
        (issue.project.id, issue.id, objectname, instrument, issue.status, ra, dec, sky, seeing, airmass, lunar_phase, transparency, estimated_hours, done_ratio))

  if level>=2 :
    # url = issue.url
    # watchers = issue.watchers
    # print('    url      =',url)
    # #print('   watchers =', end='')
    # print('    watchers =')
    # for watcher in watchers :
    #   print('              ',watcher)
      
    #for i in dir(issues[0]) :
    #  print(i)
    #  value =  issue[i]
    #if level>=3 :
    #  print('  journals')
    #  print('   ',len(issue.journals))

    if level>=3 :
      print('  project\n    %-s  %-s' % (str(issue.project.id),str(issue.project.name)))

    for dd in reversed(list(issue)) :
      keywd = dd[0]
      value = dd[1]
      #print(">>>>>>",dd,keywd)
      if level>=5 :
        print(( '  {:18s}  {:10s}').format(keywd,str(type(value))))
      if level>=2 :
        if keywd=='journals' :
          if level<=4 :
            print(' ',keywd)
          print('    len = ',len(issue.journals))
          print('    --------------------------------------------------------------------------------------------------------')
          nj=1
          for journal in issue.journals:
            if last==0 or nj<=last :
              print('   ',journal.created_on,' ',journal.user.name,'  (',journal.id,')')
              print('    --------------------------------------------------------------------------------------------------------')
              if not journal.notes is None :
                for line in journal.notes.split('\n') :
                  linenew = re.findall('.{1,100}(?:\W|$)',line)
                  for l in linenew :
                    print('      |',l)
              print('    --------------------------------------------------------------------------------------------------------')
            if not journal.notes is None :
              if len(journal.notes)>0 :
                nj+=1
        if keywd=='status' :
          if level<=4 :
            print(' ',keywd)
            print('   ',value['name'])
      if level>=3 :
        if keywd=='description' :
          if level<=4 :
            print(' ',keywd)
          for line in value.split('\n') :
            linenew = re.findall('.{1,120}(?:\W|$)',line)
            for l in linenew :
              print('    |',l)
        if keywd=='author' or keywd=='assigned_to' :
          if level<=4 :
            print(' ',keywd)
            print('   ',value['name'])
      if level>=4 :
        if keywd=='custom_fields' :
          if level<=4 :
            print(' ',keywd)
          if isinstance(value, list) :
            for j in value :
              print('    {:15s}       {}'.format(j['name'],j['value']))
          #if keywd=='journals' :
          #  if level<=4 :
          #    print(' ',keywd)
          #    print('    len = ',len(issue.journals))
      if level>=6 :
        if keywd=='time_entries' :
          print('    len = ',len(issue.time_entries))
          for time_entries in issue.time_entries:
            print('    | ',time_entries)
        elif keywd!='description' and keywd!='custom_fields' :
          if   isinstance(value, float) :
             print('   ',value)         
          elif isinstance(value, int) :
             print('   ',value)         
          elif isinstance(value, str) :
             print('   ',value)         
          elif isinstance(value,datetime.datetime) :
             print('   ',value)         
          elif isinstance(value,datetime.date) :
             print('   ',value)         
          elif isinstance(value,redminelib.managers.standard.IssueManager) :
             print('   ',value)         
          #elif isinstance(value,redminelib.resultsets.ResourceSet) :
          #  print('    = ',value)         
          elif value is None :
            print('   ',None) 
          elif isinstance(value, list) :
            for j in value :
              print('    {:15s}       {}'.format(j['name'],j['value']))
          elif isinstance(value, dict) :
            #print(value)
            print('    {:15s}       (id = {:<3})'.format(value['name'],value['id']))
  return


    
def get_project_details(redmine,projID,level,last,assigned_to,full):

  #issues = redmine.issue.all(sort='category:desc', include=['relations', 'attachments'])
  #issues = redmine.issue.all(sort='category:desc', include=['relations'])
  # issues = redmine.issue.all()
  # for issue in issues :
  #   print( issue )
  
  #projects = redmine.project.get('Observations')
  # if projID!=0 :
  project = redmine.project.get(projID)
  #  #projects = redmine.project.get(projID, include=['trackers', 'issue_categories', 'enabled_modules', 'time_entry_activities', 'issue_custom_fields'])
  # else :
  #projects = redmine.project.all()

  #iiii = redmine.issue.get(issueID)
  #print(iiii.id,iiii.estimated_hours)
  #return
  #for project in list(projects) :

  #print("=====================================")
  #print('%-6i   %-32s' % (projID, project.name))
  #print("-------------------------------------")

  print('')

  if level==0 :
    # for i in project.issues :
    #   print(i.url,list(i.watchers))
    #print(list(p.issue_custom_fields.values()))
    # print("  ------------------------------------")
    # for i in dir(project) :
    #   print(i)
    for i in list(project) :
      #`print(i)
      print( '%21s  %s' % (i[0],i[1]) )   
    print( '%21s  %s' % ('keywords',dir(project.issues[0])) )
    return
  
  # issues = redmine.issue.filter(project_id=projID, status__name='Open')
  # issues = redmine.issue.filter(project_id=projID, include=['relations', 'attachments'])
  # issues = redmine.issue.filter(project_id=projID)
  # issues = redmine.issue.filter(project_id=64)
  if full :
    issues = redmine.issue.filter(project_id=int(project.id),status_id='*')
  else : 
    issues = redmine.issue.filter(project_id=int(project.id))
    
  if level<=1 :
    print("======================================================================================================================================================")
    print('%-3s  %-7s  %-20s  %-10s  %-11s  %6s  %6s  %7s  %7s  %6s  %14s  %14s  %9s  %4s' % ( 'pID','issueID', 'Objectname', 'Instrument', 'Status', 'RA', 'DEC', 'Sky_min', 'PSF_max', 'AM_max', 'moonphase', 'transparency' ,'est.hours', 'done'))
    print("------------------------------------------------------------------------------------------------------------------------------------------------------")

  for issue in reversed(issues) :
      print_issue(issue,level,last,assigned_to)

  print("======================================================================================================================================================")
  print('')
  return
      
        
def get_issue_details(redmine,issueID,level,last,assigned_to):
  issue = redmine.issue.get(issueID)
  print('')
  if level<=1 :
    print("======================================================================================================================================================")
    print('%-3s  %-7s  %-20s  %-10s  %-11s  %6s  %6s  %7s  %7s  %6s  %14s  %14s  %9s  %4s' % ( 'pID','issueID', 'Objectname', 'Instrument', 'Status', 'RA', 'DEC', 'Sky_min', 'PSF_max', 'AM_max', 'moonphase', 'transparency' ,'est.hours', 'done'))
    print("------------------------------------------------------------------------------------------------------------------------------------------------------")
  print_issue(issue,level,last,assigned_to)
  print("======================================================================================================================================================")
  print('')
  return
      


######################################## MAIN ########################################


help = '''\

examples:
  ./redmine.py -h
  ./redmine.py
  ./redmine.py -p 37
  ./redmine.py -p 37 -a "Arno Riffeser"
  ./redmine.py -p 95
  ./redmine.py -p 95 -f
  ./redmine.py -i 3617
  ./redmine.py -i 3617 -v 2  -l 1
  ./redmine.py -i 3617 -v 2  -l 0          
  ./redmine.py -i 3617 -v 3  -l 2    # VERY USEFUL 
  ./redmine.py -i 3617 -v 4
  ./redmine.py -i 3617 -v 5
  ./redmine.py -i 3617 -v 6

important:
  use a terminal width of min. 150
'''

parser = argparse.ArgumentParser(epilog=help,formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='redmine.py  2025-06-16  (Arno Riffeser, USM)\n')
parser.add_argument('-v', dest='level',       type=int,   default=1,  help='[%(default)s] level')
parser.add_argument('-p', dest='projID',      type=int,   default=0,  help='[%(default)s] projID')
parser.add_argument('-i', dest='issueID',     type=int,   default=0,  help='[%(default)s] issueID')
parser.add_argument('-l', dest='last',        type=int,   default=2,  help='[%(default)s] nr last journals')
parser.add_argument('-a', dest='assigned_to', type=str,   default='', help='[%(default)s] assigned_to')
parser.add_argument('-f', dest='full',        action='store_true',    help='[%(default)s] full list (also closed)')
args = parser.parse_args()

APIKEY=''
if APIKEY=='' :
  print('missing API access key')
  exit(-1)
redmine = Redmine('https://luna.mpe.mpg.de/redmine',key=APIKEY)

if args.projID==0 and args.level==0 :
  get_structure(redmine)

elif args.projID==0 and args.issueID==0 :
  get_all_projects(redmine)

elif args.issueID>0 :
  get_issue_details(redmine,args.issueID,args.level,args.last,args.assigned_to)

elif args.projID>0 :
  get_project_details(redmine,args.projID,args.level,args.last,args.assigned_to,args.full)

else :
  print('wrong parameters')
