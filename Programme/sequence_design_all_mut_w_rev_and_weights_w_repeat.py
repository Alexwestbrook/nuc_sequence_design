#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 14:22:48 2019

@author: epierre
"""
import numpy as np
from math import exp, sqrt
import matplotlib.pyplot as plt
import os
import argparse

try:
    from keras.models import load_model
except ModuleNotFoundError:
    from tensorflow.keras.models import load_model
import gc

from MyModuleLibrary.array_modifier import rolling_window
from MyModuleLibrary.mykeras.losses import correlate, mse_var, mae_cor


class Sequence:

    def __init__(self, length, repeat):
        """
            Initiate a DNA sequence

            Args:
                length: the length of the artificial repeated element.
                repeat: the number of repetition of the sequence on which to
                predict (1 is a good number...)
        """
        self.length = length
        self.repeat = repeat
        '''pure random sequence 25% each'''
        '''sequence with weights on singlets'''
        '''yeast-like sequence'''
        self.sequence = self._get_yeast_like_sequence()
        '''pre-saved sequence'''
        '''reverse complementary of a pre-saved sequence'''
        self.mutated_seq = self._duplicate(self.sequence)
        self.atgc_content = self._get_atgc_content(self.mutated_seq)
        self.GC_content = self._get_gc_content_all(self.mutated_seq)
        self.GC_energy = self._get_gc_energy(self.GC_content)
        self.WX = 2001

    @property
    def seq_predictable(self):
        """
            Change the sequence into a serie of one-hot-encoded window to serve
            as input of a model.
        """
        x_seq = self._process(self.sequence)
        return x_seq

    @property
    def mutated_seq_predictable(self):
        """
            Change one nucleotid in the mutated_sequence and convert it into a
            predictable form.
        """
        x_mutated_seq = self._process_mutated(self.mutated_seq)
        return x_mutated_seq

    def _one_hot_encoder(self, nucleotid):
        res = (np.arange(nucleotid.max()) == nucleotid[..., None]-1
               ).astype(int)
        res = res.reshape(res.shape[0], 4)
        return res

    def _process(self, nucleotid):
        x = self._duplicate(nucleotid)
        self.mutated_seq = self._all_mutation(x)
        self.GC_content = self._get_gc_content_all(self.mutated_seq)
        self.GC_energy = self._get_gc_energy(self.GC_content)
        x_w = self._get_reverse_complementary(self.mutated_seq)
        x_r = self._rescale(x_w)
        x_r = x_r.flatten()
        x_r = self._one_hot_encoder(x_r)
        HALF_WX = self.WX // 2
        nw_seq_length = self.repeat * self.length + 2 * HALF_WX
        x_r = x_r.reshape(6*self.length, nw_seq_length, 4)
        x_slide = rolling_window(x_r, window=(self.WX, 4))
        x_seq = x_slide.reshape(x_slide.shape[0], x_slide.shape[1],
                                self.WX, 4, 1)
        x_seq = x_seq.reshape(6 * self.length * self.length * self.repeat,
                              self.WX, 4, 1)
        return x_seq

    def _process_mutated(self, nucleotid):
        self.mutated_seq = self._all_mutation(nucleotid)
        self.GC_content = self._get_gc_content_all(self.mutated_seq)
        self.GC_energy = self._get_gc_energy(self.GC_content)
        x_w = self._get_reverse_complementary(self.mutated_seq)
        x_r = self._rescale(x_w)
        x_r = x_r.flatten()
        x_r = self._one_hot_encoder(x_r)
        HALF_WX = self.WX // 2
        nw_seq_length = self.repeat * self.length + 2 * HALF_WX
        x_r = x_r.reshape(6*self.length, nw_seq_length, 4)
        x_slide = rolling_window(x_r, window=(self.WX, 4))
        x_seq = x_slide.reshape(x_slide.shape[0], x_slide.shape[1],
                                self.WX, 4, 1)
        x_seq = x_seq.reshape(6*self.length*self.length*self.repeat,
                              self.WX, 4, 1)
        return x_seq

    def _rescale(self, nucleotid):
        HALF_WX = self.WX // 2
        MARGIN = 2 * HALF_WX
        repeat_number = MARGIN // self.length + self.repeat + 2
        start = (HALF_WX // self.length + 1) * self.length - HALF_WX
        stop = start + self.repeat * self.length + 2 * HALF_WX
        new_sequence = np.tile(nucleotid, (1, repeat_number))
        return new_sequence[:, start:stop]

    def _all_mutation(self, nucleotid):
        """mute 1 bp per line in a circular rotation 1 -> 2 -> 3 -> 4 -> 1"""
        for i in range(3*self.length):
            nucleotid[i, i//3] = (nucleotid[i, i//3]
                                  + 1
                                  + i % 3
                                  - 4 * ((nucleotid[i, i//3] + i % 3) // 4))
            # Equivalent to (nucleotid[i, i//3] + i % 3) % 4 + 1
        return nucleotid

    def _duplicate(self, nucleotid):
        x_lin = np.ravel(nucleotid)
        x_dup = np.tile(x_lin, (3*(self.length), 1))
        return x_dup

    def _get_atgc_content(self, nucleotid):
        counter = np.zeros((4, 1))
        for i in range(self.length):
            counter[nucleotid[0, i]-1, 0] += 1
        counter /= self.length
        return counter

    def _get_gc_content_all(self, nucleotid):
        counter = np.zeros((3*(self.length), 1))
        for i in range(3*self.length):
            for j in range(self.length):
                if nucleotid[i, j] > 2.5:
                    counter[i, 0] += 1
        counter /= self.length
        return counter

    def _get_gc_energy(self, gc_content):
        coeff = 1.
        natural = 0.3806
        GC_energy = coeff*np.sqrt((gc_content-natural)*(gc_content-natural))
        GC_energy_1D = GC_energy.reshape((3*self.length,))
        return GC_energy_1D

    def _get_atgc_content_mutated(self, nucleotid):
        counter = np.zeros((4, 1))
        for i in range(self.length):
            counter[nucleotid[i]-1, 0] += 1
        counter /= self.length
        return counter

    def _get_reverse_complementary(self, nucleotid):
        complementary_nucleotid = np.copy(nucleotid)
        for i in range(len(complementary_nucleotid)):
            complementary_nucleotid[i] = (complementary_nucleotid[i]
                                          + 1
                                          - 2*(complementary_nucleotid[i]//2)
                                          + 2*(complementary_nucleotid[i]//3))
        reverse_complementary_nucleotid = np.flip(complementary_nucleotid, 1)
        whole_seq = np.concatenate((nucleotid,
                                    reverse_complementary_nucleotid))
        return whole_seq

    def _get_pure_random_sequence(self):
        L = self.length
        seq = np.random.randint(1, 5, (L, 1))
        return seq

    def _get_weigth_singlet_sequence(self):
        choices = [1, 2, 3, 4]  # A, T, G, C
        weights = [0.310052253746, 0.309303360736,
                   0.19039708206, 0.190247303458]
        L = self.length
        seq = np.random.choice(choices, p=weights, size=(L, 1))
        return seq

    def _get_yeast_like_sequence(self):
        dconv = {'A': 1, 'T': 2, 'G': 3, 'C': 4}
        items_triplet = [["AAA", 0.0273654778875],
                         ["AAT", 0.0307494815188],
                         ["AAG", 0.0228417497464],
                         ["AAC", 0.0188620559663],
                         ["ATA", 0.0232890668877],
                         ["ATT", 0.0306549624177],
                         ["ATG", 0.0190686791175],
                         ["ATC", 0.0183432999696],
                         ["AGA", 0.0202446725847],
                         ["AGT", 0.0157385293927],
                         ["AGG", 0.0120885534071],
                         ["AGC", 0.011854453773],
                         ["ACA", 0.0171112544773],
                         ["ACT", 0.0159011901713],
                         ["ACG", 0.0088320406564],
                         ["ACC", 0.011652226859],
                         ["TAA", 0.0232154299135],
                         ["TAT", 0.0232176280322],
                         ["TAG", 0.0135261229913],
                         ["TAC", 0.0147405835346],
                         ["TTA", 0.0233297320823],
                         ["TTT", 0.0272434823035],
                         ["TTG", 0.0241617199839],
                         ["TTC", 0.0247321317684],
                         ["TGA", 0.0210414905882],
                         ["TGT", 0.0169826645374],
                         ["TGG", 0.015275825421],
                         ["TGC", 0.0130623199603],
                         ["TCA", 0.0209546649023],
                         ["TCT", 0.0198698933583],
                         ["TCG", 0.00931562675506],
                         ["TCC", 0.0132887261793],
                         ["GAA", 0.0249464483349],
                         ["GAT", 0.0183993519946],
                         ["GAG", 0.011209305955],
                         ["GAC", 0.010041005903],
                         ["GTA", 0.0145614368662],
                         ["GTT", 0.0188235888903],
                         ["GTG", 0.0103828133501],
                         ["GTC", 0.00993000091222],
                         ["GGA", 0.0132788346454],
                         ["GGT", 0.0118423641205],
                         ["GGG", 0.00563377804937],
                         ["GGC", 0.00801653864457],
                         ["GCA", 0.0128513005718],
                         ["GCT", 0.0120588788056],
                         ["GCG", 0.00559531097334],
                         ["GCC", 0.00818469471979],
                         ["CAA", 0.0242914089831],
                         ["CAT", 0.0189851506096],
                         ["CAG", 0.0130161594691],
                         ["CAC", 0.0102003695037],
                         ["CTA", 0.0135239248727],
                         ["CTT", 0.0227450325267],
                         ["CTG", 0.0130535274858],
                         ["CTC", 0.011090607549],
                         ["CGA", 0.00936398536493],
                         ["CGT", 0.00882984253777],
                         ["CGG", 0.00577335858239],
                         ["CGC", 0.00575577363335],
                         ["CCA", 0.0152285658705],
                         ["CCT", 0.0119160010947],
                         ["CCG", 0.00598108079295],
                         ["CCC", 0.0059338212424]]

        elems_triplet = [i[0] for i in items_triplet]
        probs_triplet = [i[1] for i in items_triplet]
        trials = self.length//3
        res = np.random.choice(elems_triplet, trials, p=probs_triplet)
        conc_res = ''.join(map(str, res))
        i = 0
        L = self.length
        seq = np.ones((L, 1), dtype=int)
        while i < 3*trials:
            seq[i] = int(dconv[conc_res[i]])
            i += 1
        return seq

    def _load_sequence(self):
        seq = np.load(os.path.join(os.path.dirname(__file__),
                                   '..',
                                   'Results_nucleosome',
                                   'input_test_seq.npy'))
        return seq

    def _rev_compl_loaded_seq(self):
        ini_seq = np.load(os.path.join(os.path.dirname(__file__),
                                       '..',
                                       'Results_nucleosome',
                                       'input_test_seq.npy'))
        complementary_nucleotid = np.copy(ini_seq)
        for i in range(len(complementary_nucleotid)):
            complementary_nucleotid[i] = (complementary_nucleotid[i]
                                          + 1
                                          - 2*(complementary_nucleotid[i]//2)
                                          + 2*(complementary_nucleotid[i]//3))
        seq = np.flip(complementary_nucleotid, 0)
        return seq


class Energy:

    def __init__(self, seq, model, y_target, length, repeat):
        """
            This class is aimed at storing the previous energy so that to be
            able to reject a mutation.

            To initialise the energy we need the model on which the prediction
            will be made, the target function and also the first sequence (in a
            predictable shape). For representation purpose we also keep the
            predicted nucleosome density.

            Args:
                seq: numpy array of the one-hot-encoded initial sequence
                (with a rolling window applied so that prediction can be made).
                model: the trained keras model (could be multi-output)
                y_target: the target nucleosome density.
            Returns:
                energy: an Energy instance.
        """
        self.model = model
        self.y_target = y_target
        self.length = length
        self.repeat = repeat
        self.reverse_y_target = self._get_reverse_target(self.y_target)
        self.prediction = self.model.predict(seq)

        self.prediction_profile = self._get_prediction_profile(self.prediction)
        (self.energy_profile,
         self.reg_energy_profile,
         self.rev_energy_profile
         ) = self._get_energy_profile(self.prediction_profile)

    def compute_mutated_energy(self, mutated_seq):
        """
            Calculate the energy of a mutated sequence and returns it.
            Args:
                mutated_seq: the one-hot-encoded proposed mutated sequence.
            Returns:
                Change self.mutated_energy to the new value.
        """
        self.mutated_pred = self.model.predict(mutated_seq)
        self.prediction_profile = self._get_prediction_profile(
            self.mutated_pred)
        (self.energy_profile,
         self.reg_energy_profile,
         self.rev_energy_profile
         ) = self._get_energy_profile(self.prediction_profile)

    def _get_corr(self, y):
        X = y - np.mean(y)
        Y = self.y_target - np.mean(self.y_target)
        sigma_XY = np.sum(X*Y)
        sigma_X = np.sqrt(np.sum(X*X))
        sigma_Y = np.sqrt(np.sum(Y*Y))
        return sigma_XY/(sigma_X*sigma_Y + (7./3 - 4./3 - 1))  # () is epsilon

    def _get_reverse_corr(self, y):
        X = y - np.mean(y)
        Y = self.reverse_y_target - np.mean(self.reverse_y_target)
        sigma_XY = np.sum(X*Y)
        sigma_X = np.sqrt(np.sum(X*X))
        sigma_Y = np.sqrt(np.sum(Y*Y))
        return sigma_XY/(sigma_X*sigma_Y + (7./3 - 4./3 - 1))  # () is epsilon

    def _get_one_energy(self, y):
        return 5*np.mean(np.abs(y - self.y_target)) - self._get_corr(y) + 1

    def _get_one_reverse_energy(self, y):
        return (5*np.mean(np.abs(y - self.reverse_y_target))
                - self._get_reverse_corr(y) + 1)

    def _get_energy(self, y):
        return np.mean([self._get_one_energy(y[:, i])
                        for i in range(y.shape[1])])

    def _get_reverse_energy(self, y):
        return np.mean([self._get_one_reverse_energy(y[:, i])
                        for i in range(y.shape[1])])

    def _get_prediction_profile(self, y):
        y1 = int(6*self.length)
        y2 = int(self.length*self.repeat)
        return y.reshape(y1, y2, 1)

    def _get_energy_profile(self, y):
        mid_length = y.shape[0]//2
        reg_energy = np.array([self._get_energy(y[i, :])
                               for i in range(mid_length)])
        rev_energy = np.array([self._get_reverse_energy(y[i, :])
                               for i in range(mid_length, y.shape[0])])
        full_energy = reg_energy + rev_energy
        return full_energy, reg_energy, rev_energy

    def _get_reverse_target(self, y):
        reverse_y = np.flip(y, 0)
        return reverse_y


def step(sequence, energy, temp, mut_profile, weight):
    """
        Propose a mutation of the sequence and accept the best mutation
        on the basis of Kinetic Monte Carlo algorithm.

        We propose a mutation and calculate the predicted density of nucleosome
        on the mutated sequence. The energy is the euclidian distance between
        the density predict and the density that we finally want to have.

        Args:
            sequence: a Sequence instance
            energy: Energy instances corresponding to the sequence
            temp: the temperature of the simulation
            mut_profile: contain the list of already done mutations. Used as
                         weight
            weight: user-defined weight for GC, reg_seq, rev_seq and
                    mut_profile
        Return:
            index: the index of the selected mutation
    """
    GC_energy = sequence.GC_energy
    reg_energy = energy.reg_energy_profile
    rev_energy = energy.rev_energy_profile
    mutation_energy = mut_profile

    weight_GC = float(weight[0])
    weight_reg = float(weight[1])
    weight_rev = float(weight[2])
    weight_mut = float(weight[3])

    used_energy = (weight_GC * GC_energy
                   + weight_reg * reg_energy
                   + weight_rev * rev_energy
                   + weight_mut * mutation_energy)
    esum = np.sum(np.exp(-(used_energy)/temp))
    print('esum:', esum)
    rand = np.random.uniform(0, esum)
    print('rand:', rand)
    index = 0
    eless = 0
    while eless <= rand:
        eless = eless + exp(-(used_energy[index])/temp)
        index = index + 1
    print('chosen index:', index-1)
    return index-1


def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--length', default=167,
                        help='''the length of the artificial sequence''')
    parser.add_argument('-r', '--repeat', default=1,
                        help='''the number of artificial sequence we want to
                        predict on''')
    parser.add_argument('-s', '--steps', default=10000,
                        help='''Number of KMC iterations.''')
    parser.add_argument('-m', '--model',
                        default="weights_with_rev_compl_rep2.hdf5",
                        help='''Name of the model to be used for prediction''')
    parser.add_argument('-t', '--temperature', default=0.1,
                        help='''temperature used in the KMC algorithm''')
    parser.add_argument('-d', '--directory', default="test",
                        help='''directory where the optput sequences will be
                        saved''')
    parser.add_argument('-w', '--weight', nargs=4, default=[1, 1, 1, 1],
                        help='''set of 4 weights used for optimization''')
    parser.add_argument('-g', '--gauss', nargs=2, default=[0.4, 0.2],
                        help='''amplitude and background of the gaussian
                                target''')
    return parser.parse_args(args)


def make_output_folder(folder):
    path_to_new_folder = os.path.join('..', 'Results_nucleosome', folder)
    visu_fold = folder + '_plot'
    path_to_visu = os.path.join('..', 'Results_nucleosome', visu_fold)
    try:
        # Create target Directory
        os.mkdir(path_to_new_folder)
        print("Directory ", folder,  " Created ")
        os.mkdir(path_to_visu)
        print("Directory ", visu_fold,  " Created ")
    except OSError:
        print("Directory ", folder,  " already exists")
        print("Directory ", visu_fold,  " already exists")


def main(command_line_arguments=None):
    """
        KMC algorithm to design a sequence that leads to a desired
        nucleosome density.
    """

    np.random.seed(2)  # for debugg

    args = _parse_arguments(command_line_arguments)

    model = load_model(
        os.path.join(
            os.path.dirname(__file__),
            '..',
            'Results_nucleosome',
            os.path.basename(os.path.abspath(args.model))),
        custom_objects={'correlate': correlate,
                        'mse_var': mse_var,
                        'mae_cor': mae_cor})

    amp = float(args.gauss[0])  # 0.4 # 0.3
    shift = float(args.gauss[1])  # 0.2 # 0.2
    len601 = 147
    x_gauss = np.arange(1, len601, 1)
    y_target_ = shift + amp*np.exp(-((x_gauss - ((len601-1)/2.))**2
                                     / (2.*len601*len601/16.)))
    # linkers at last value
    y_target_ = np.append(y_target_,
                          np.repeat(y_target_[-1],
                                    (int(args.length) - (len601-1))))
    y_target = y_target_
    for i in range(int(args.repeat) - 1):
        y_target = np.append(y_target, y_target_)

    make_output_folder(args.directory)

    sequence = Sequence(int(args.length), int(args.repeat))
    energy = Energy(sequence.seq_predictable, model, y_target,
                    int(args.length), int(args.repeat))
    print('Initialization done!')
    counter = 0

    plt.ion()
    fig = plt.figure()
    ax_pred = fig.add_subplot(1, 1, 1)
    plt.xlabel('bp')
    plt.ylabel('Nucleosome occupancy')
    line_pred, = ax_pred.plot(y_target, 'r', label='target')
    mean_prediction = np.mean(energy.prediction_profile[0, :, :], axis=1)
    mean_com_prediction = np.mean(
        energy.prediction_profile[3*sequence.length, :, :], axis=1)
    mean_rev_prediction = np.flip(mean_com_prediction, 0)  # flipped
    line_pred_reg, = ax_pred.plot(mean_prediction, 'b', label='pred')
    line_pred_rev, = ax_pred.plot(mean_rev_prediction, 'k',
                                  label='rev_com_pred')
    ax_pred.legend()
    count_mutation = np.zeros((int(3*sequence.length),))
    energy_file = open(os.path.join(os.path.dirname(__file__),
                                    '..',
                                    'Results_nucleosome', args.directory,
                                    'energy.txt'), 'w')

    with open(os.path.join(os.path.dirname(__file__),
                           '..',
                           'Results_nucleosome', args.directory,
                           'config.txt'), 'w') as f_head:
        f_head.write("length\trepeat\tstep\tmodel\ttemperature\tdirectory"
                     "\tweight_GC\tweight_seq\tweight_rev_comp\tweight_mut"
                     "\ttarget_amplitude\ttarget_background\n")
        len_w = str(args.length)
        rep_w = str(args.repeat)
        stp_w = str(args.steps)
        mod_w = str(args.model)
        temp_w = str(args.temperature)
        dir_w = str(args.directory)
        GC_w = str(args.weight[0])
        seq_w = str(args.weight[1])
        rev_w = str(args.weight[2])
        mut_w = str(args.weight[3])
        amp_w = str(args.gauss[0])
        back_w = str(args.gauss[1])
        f_head.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'
                     % (len_w, rep_w, stp_w, mod_w, temp_w, dir_w, GC_w,
                        seq_w, rev_w, mut_w, amp_w, back_w))

    dir_plot = str(args.directory)+'_plot'

    for i in range(int(args.steps)):
        accepted = step(sequence, energy, float(args.temperature),
                        count_mutation, args.weight)
        sequence.sequence = sequence.mutated_seq[accepted, :]
        mean_prediction = np.mean(energy.prediction_profile[accepted, :, :],
                                  axis=1)
        mean_com_prediction = np.mean(
            energy.prediction_profile[accepted + 3*sequence.length, :, :],
            axis=1)
        mean_rev_prediction = np.flip(mean_com_prediction, 0)  # flipped
        count_mutation[accepted] = count_mutation[accepted] + 1
        reg_nrj = energy.reg_energy_profile[accepted]
        rev_nrj = energy.rev_energy_profile[accepted]
        GC_nrj = sequence.GC_energy[accepted]
        full_nrj = reg_nrj + rev_nrj + GC_nrj
        s = str(i) + ', ' + str(full_nrj)+'\n'
        energy_file.write(s)

        if accepted:
            counter += 1
            print('nb of iteration:', counter)

        np.save(os.path.join(os.path.dirname(__file__),
                             '..',
                             'Results_nucleosome', args.directory,
                             'designed_sequence_'+str(i)+'.npy'),
                sequence.sequence)

        if i % 1 == 0:  # plot every 1 step
            line_pred_reg.set_ydata(mean_prediction)
            line_pred_rev.set_ydata(mean_rev_prediction)
            ax_pred.set_ylim(min(np.min(mean_prediction),
                                 np.min(y_target)) - 0.05,
                             max(np.max(mean_prediction),
                                 np.max(y_target)) + 0.05)

            fig.canvas.draw()
            plt.savefig(os.path.join(os.path.dirname(__file__),
                                     '..',
                                     'Results_nucleosome', dir_plot,
                                     'plot_step_'+str(i)+'.png'))
            plt.pause(1e-17)

        if all(i > 0 for i in count_mutation):
            count_mutation = np.zeros((int(3*sequence.length),))
        energy.compute_mutated_energy(sequence.seq_predictable)
        _ = gc.collect()  # fixes memory leak from model.predict
    plt.ioff()
    plt.show()

    energy_file.close()


if __name__ == '__main__':
    main()
