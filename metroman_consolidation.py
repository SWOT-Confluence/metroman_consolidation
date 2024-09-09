#!/usr/bin/env python3
"""
Author : Travis Simmons
Date   : July 31, 2024
Purpose: Rock the Casbah
"""
# Sample deployment

# python3 prepare.py -i /home/u24/travissimmons/cjx/season10/50_hand_label_test_2020_03_02 -o /home/u24/travissimmons/cjx/season10/gifs
# makeflow process.makeflow -j 1
import argparse
import os
import glob
import json
import netCDF4 as ncf
import numpy as np


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--index',
                        help='index of continent to process',
                        metavar='int',
                        type=int,
                        )


    return parser.parse_args()


class process_reaches():

    def __init__(self, indir, input_mnt_path, index, output_dir):
        self.indir = indir
        self.input_mnt_path = input_mnt_path
        self.index = index
        self.output_dir = output_dir
        self.cont, self.cont_number = self.get_cont_info()
        self.all_reach_ids = self.find_reach_ids()
        self.parse_data()

    def get_cont_info(self):
        with open(os.path.join(self.input_mnt_path,'continent.json')) as f:
            cont_data = json.load(f)
        cont = list(cont_data[self.index].keys())[0]
        cont_number = cont_data[self.index][cont]
        return cont, cont_number

    def find_reach_ids(self):
        all_output_files = glob.glob(os.path.join(self.indir, f'{self.cont_number}*.nc'))
        all_reach_ids = []
        all_reach_ids = [os.path.basename(i).split('_')[0].split('-') for i in all_output_files]
        all_reach_ids = list(set(sum(all_reach_ids, [])))
        return all_reach_ids

    def find_reach_nt(self, reach_id):
                # open swot file for this reach
        swotfname=os.path.join(self.input_mnt_path, 'swot', f'{reach_id}_SWOT.nc')
        swotfile=ncf.Dataset(swotfname)
        
        # get swot file times
        self.tswot=list(np.round(swotfile['reach']['time'][:].filled(np.nan)/3600.))
        ntswot=len(self.tswot)    
        print('   there are ',ntswot,'times in the swot file')
        
        tswot_str=list(swotfile['reach']['time_str'][:].filled(np.nan))
        
        swotfile.close()
        
        # # set up the Q estimates
        # Q=np.full( (ntswot,),np.nan )
        # qu=np.full( (ntswot,),np.nan )
        # # print('tswot', tswot, tsf)
        # for i,t in enumerate(tswot):
        #     if not np.isnan(t) and t in tsf:
        #         idx=tsf.index(t)
        #         Q[i]=Qsetfile[ii,idx]
        #         qu[i]=qu_setfile[ii,idx]
    
    def generate_outfile(self,reach_id, data_dict ):
        # output     
        outfile= os.path.join(self.output_dir, f'{reach_id}_metroman.nc')     
        dsout = ncf.Dataset(outfile, 'w', format="NETCDF4") #output dataset
        # dsout.createDimension("nr", nr)
        dsout.createDimension("nt", self.reach_nt)

        fillvalue = np.nan
        
        t = dsout.createVariable("nt","f8",("nt"),fill_value=fillvalue)
        t.long_name= 'swot timeseries "time" variable converted to hours and rounded to integer'
        t[:] = self.tswot

        for key in list(data_dict[reach_id].keys()):
            set_group = dsout.createGroup(key)

            #                 'A0hat': A0hat,
            #     'nahat': nahat,
            #     'x1hat': x1hat,
            #     'qu_setfile': qu_setfile,
            #     'Qsetfile':Qsetfile
            # }

            
            allq = set_group.createVariable("allq", "f8", ( "nt"), fill_value=fillvalue)
            allq[:] = data_dict[reach_id][key]['Qsetfile']
            A0 = set_group.createVariable("A0hat", "f8", fill_value=fillvalue)
            A0[:]= data_dict[reach_id][key]['A0hat']
            na = set_group.createVariable("nahat", "f8", fill_value=fillvalue)
            na[:]=data_dict[reach_id][key]['nahat']
            x1 = set_group.createVariable("x1hat", "f8",  fill_value=fillvalue)
            x1[:]=data_dict[reach_id][key]['x1hat']
            allqu = set_group.createVariable("q_u", "f8", ("nt"), fill_value=fillvalue)
            allqu[:]=data_dict[reach_id][key]['qu_setfile']

        return dsout 

    def extract_mm_data(self, a_file, reach_id):
        mm_data=ncf.Dataset(a_file)
        present_reaches = os.path.basename(a_file).split('_')[0].split('-')
        reach_index = os.path.basename(a_file).split('_')[0].split('-').index(str(reach_id))

        # nr=np.max(mm_data['nr'][:].data)
        nr=mm_data.dimensions['nr'].size

        # determine number of times in the set file
        tsf=list(mm_data['t'][:].data) #setfile number of times - same for all reaches    
        ntsf=len(tsf)



        # grab the discharge array
        Qsetfile=mm_data['allq'][reach_index][:].filled(np.nan)
        Qsetfile

        # read other key metroman data
        A0hat=mm_data['A0hat'][reach_index].filled(np.nan)
        nahat=mm_data['nahat'][reach_index].filled(np.nan)
        x1hat=mm_data['x1hat'][reach_index].filled(np.nan)
        qu_setfile=mm_data['q_u'][reach_index][:].filled(np.nan)

        print(Qsetfile)
        
        mm_data.close()
        self.data_dict[reach_id][os.path.basename(a_file).split('_')[0]] = {
                'A0hat': A0hat,
                'nahat': nahat,
                'x1hat': x1hat,
                'qu_setfile': qu_setfile,
                'Qsetfile':Qsetfile
            }
        # self.data_dict[reach_id]['average'] = {
        #         'A0hat': self.data_dict[reach_id]['average']['A0hat'].append([A0hat]),
        #         'nahat': self.data_dict[reach_id]['average']['nahat'].append([nahat]),
        #         'x1hat': self.data_dict[reach_id]['average']['x1hat'].append([x1hat]),
        #         'qu_setfile': self.data_dict[reach_id]['average']['qu_setfile'].append(qu_setfile),
        #         'Qsetfile':self.data_dict[reach_id]['average']['Qsetfile'].append(Qsetfile)
        #     }
        # return A0hat, nahat, x1hat, Qsetfile, qu_setfile

    def create_average_group(self):
        for a_set in list(self.data_dict[self.reach_id].keys()):
            if a_set != 'average':
                for a_var in list(self.data_dict[self.reach_id][a_set].keys()):
                    self.data_dict[self.reach_id]['average'][a_var].append(self.data_dict[self.reach_id][a_set][a_var])
        for a_var in list(self.data_dict[self.reach_id]['average'].keys()):
            if len(self.data_dict[self.reach_id]['average'][a_var])>1:
                if not np.isnan(self.data_dict[self.reach_id]['average'][a_var]).all():
                    average_array_axis0 = np.nanmean(self.data_dict[self.reach_id]['average'][a_var][:], axis=0)
                    self.data_dict[self.reach_id]['average'][a_var] = average_array_axis0
                else:
                    self.data_dict[self.reach_id]['average'][a_var] = self.data_dict[self.reach_id]['average'][a_var][0]
            else:
                self.data_dict[self.reach_id]['average'][a_var] = self.data_dict[self.reach_id]['average'][a_var][0]

    def parse_data(self):


        for reach_id in self.all_reach_ids:
            self.reach_id = reach_id

            self.data_dict = {

            }
            self.data_dict[reach_id] = {}
            self.data_dict[reach_id]['average'] = {
                'A0hat': [],
                'nahat': [],
                'x1hat': [],
                'qu_setfile': [],
                'Qsetfile': []
            } 
            self.reach_nt = self.find_reach_nt(reach_id = reach_id)

            
            self.reach_outfiles = glob.glob(os.path.join(self.indir, f'*{reach_id}*'))

            for a_file in self.reach_outfiles:
                self.extract_mm_data(a_file = a_file, reach_id = reach_id)
            
            self.create_average_group()

            outfile = self.generate_outfile(reach_id = reach_id, data_dict = self.data_dict)
            
                            # run mike's thing to filter the output into the correct days
            # create_data_group_and_add_data(outfile = outfile)

            # create_average_group()

            outfile.close()


def main():
    """Make a jazz noise here"""
    args = get_args() 
    indir = '/mnt/data/flpe'
    input_mnt_path = '/mnt/data/input'
    output_dir = '/mnt/data/flpe'
    index = args.index
    process_reaches(indir, input_mnt_path, index, output_dir)

# --------------------------------------------------
if __name__ == '__main__':
    main()