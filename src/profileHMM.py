from math import ceil
import argparse
import os


def _get_states(i):  # get state from parameter
    return 'M{}'.format(i), 'I{}'.format(i), 'D{}'.format(i)


def _keyify(i):  # for sort
    return int("{}".format(i[0][1:]))


class ProfileHiddenMarkovMoldel:

    def __init__(self, inputfile, output):
        self.inputfile = inputfile
        self.output = output
        self.t_prob = {}  # transition probability
        self.e_prob = {}  # emission probability
        self.inp_strings = []  # input string
        self.char_list = set()

        # Open file and read input strings
        with open(inputfile, 'r') as f:
            for line in f.read().splitlines():
                self.inp_strings.append(line)
                self.char_list = self.char_list.union(set(line))
        self.char_list = self.char_list - {'.'}  # Clean the data

        # number of input strings
        self.num_of_strings = len(self.inp_strings)
        # number of char in a string
        self.num_of_chars = len(self.inp_strings[0])

        self.frequncy_list = [{} for i in range(self.num_of_chars + 1)]
        for string in self.inp_strings:
            for index, char in enumerate(string):
                if char in self.frequncy_list[index].keys():
                    self.frequncy_list[index][char] += 1
                else:
                    self.frequncy_list[index][char] = 1

        # Which states are the match state
        self.match_states = [
            k for n, k in zip(self.frequncy_list, range(self.num_of_chars + 1))
            if int(n.get('.', 0)) < ceil(self.num_of_strings / 2)
        ]

        # State lists(temp)
        match_state = ['M{}'.format(k)
                       for k in range(0, len(self.match_states) + 1)]
        insert_state = ['I{}'.format(k)
                        for k in range(0, len(self.match_states))]
        delete_state = ['D{}'.format(k)
                        for k in range(1, len(self.match_states))]

        # formatting the transition probabilities
        self.t_prob.update({key: {'strs': []} for key in match_state})
        self.t_prob.update({key: {'strs': []} for key in insert_state})
        self.t_prob.update({key: {'strs': []} for key in delete_state})

        # put all input chars(proteins) at the begining state
        self.t_prob['M0']['strs'] = [n for n in range(self.num_of_strings)]

    def build_model(self):
        i = 0  # counter
        j = 0  # current state no
        while i < self.num_of_chars + 1:
            M, I, D, = _get_states(j)
            nextM, nextD = _get_states(j + 1)[::2]

            # If the current state is the match state
            if i in self.match_states:
                # D --> D list and D --> M list
                deltodel, deltomatch = [], []
                # I --> M List and I --> D List
                instomatch, instodel = [], []
                # M --> D List and M --> M List
                matchtodel, matchtomatch = [], []

                # D --> D and D --> M
                # D(j) --> D(j+1) or D(j) --> M(j+1)
                if self.t_prob.get(D, {}).get('strs', []) and i != 0:
                    try:
                        deltodel = [n for n in self.t_prob[D]['strs']
                                    if self.inp_strings[n][i] == '.']
                    except IndexError:
                        pass
                    deltomatch = [
                        n for n in self.t_prob[D]['strs'] if n not in deltodel
                    ]

                    # If deltodel is not empty
                    # D --> D
                    if deltodel:
                        self.t_prob[D][nextD] = {
                            'prob': float(len(deltodel) /
                                          len(self.t_prob[D]['strs'])),
                            'strs': deltodel
                        }
                        self.t_prob[nextD]['strs'].extend(deltodel)

                    # If deltomatch is not empty
                    # D --> M
                    if deltomatch:
                        self.t_prob[D][nextM] = {
                            'prob': float(len(deltomatch) /
                                          len(self.t_prob[D]['strs'])),
                            'strs': self.t_prob[D]['strs']
                        }
                        self.t_prob[nextM]['strs'].extend(deltomatch)

                # I --> M and I --> D
                # I(j) --> D(j+1) or I(j) --> M(j+1)
                if self.t_prob[I]['strs'] and i != 0:
                    try:
                        instodel = list(set([n for n in self.t_prob[I]['strs']
                                             if self.inp_strings[n][i] == '.']))
                    except IndexError:
                        pass
                    instomatch = list(set(
                        [n for n in self.t_prob[I]['strs'] if n not in instodel]))

                    # if instodel is not empty
                    # I --> D
                    if instodel:
                        self.t_prob[I][nextD] = {
                            'prob': float(len(instodel) /
                                          len(self.t_prob[I]['strs'])),
                            'strs': instodel
                        }
                        self.t_prob[nextD]['strs'].extend(set(instodel))

                    # If instomatch is not empty
                    # I --> M
                    if instomatch:
                        self.t_prob[I][nextM] = {
                            'prob': float(len(instomatch) /
                                          len(self.t_prob[I]['strs'])),
                            'strs': instomatch
                        }
                        self.t_prob[nextM]['strs'].extend(set(instomatch))

                # M --> D and M --> M
                # M(j) --> D(j+1) or M(j) --> M(j+1)
                if self.t_prob[M]['strs']:
                    try:
                        matchtodel = [n for n in self.t_prob[M]['strs']
                                      if self.inp_strings[n][i] == '.'
                                      and n not in self.t_prob[I]['strs']]
                    except IndexError:
                        pass

                    matchtomatch = [
                        n for n in self.t_prob[M]['strs']
                        if n not in matchtodel + self.t_prob[I]['strs']
                    ]

                    # If matchtodel is not empty
                    # M --> D
                    if matchtodel:
                        self.t_prob[M][nextD] = {
                            'prob': float(len(matchtodel) /
                                          len(self.t_prob[M]['strs'])),
                            'strs': matchtodel
                        }
                        self.t_prob[nextD]['strs'].extend(matchtodel)

                    # If matchtomatch is not empty
                    # M --> M
                    if matchtomatch:
                        self.t_prob[M][nextM] = {
                            'prob': float(len(matchtomatch) /
                                          len(self.t_prob[M]['strs'])),
                            'strs': matchtomatch
                        }
                        self.t_prob[nextM]['strs'].extend(matchtomatch)
                j += 1
            else:
                insert_states = []

                # This loop for going to the next insert state
                while True:
                    insert_states.extend([n for n in range(self.num_of_strings)
                                          if self.inp_strings[n][i] != '.'])
                    if i + 1 in self.match_states or i + 1 == self.num_of_chars:
                        # if i+1 is not match state or i+1 last char in strings
                        break
                    i += 1  # next insert state

                # If the current insert_state is not empty
                if insert_states:
                    # M --> I
                    come_from_match = [n for n in self.t_prob[M]['strs']
                                       if n in insert_states]
                    # D --> I
                    come_from_del = [n for n in self.t_prob.get(D, {})
                                     .get('strs', []) if n in insert_states]
                    # I --> I
                    come_from_ins = [n for n in set(insert_states) for k in
                                     range(insert_states.count(n) - 1)]

                    # If the string comes from the match state
                    # M(j) --> I(j)
                    if come_from_match:
                        self.t_prob[M][I] = {
                            'prob': float(len(come_from_match) /
                                          len(self.t_prob[M]['strs'])),
                            'strs': come_from_match
                        }

                    # If the string comes from the delete state
                    # D(j) --> I(j)
                    if come_from_del:
                        self.t_prob[D][I] = {
                            'prob': float(len(come_from_del) /
                                          len(self.t_prob[D]['strs'])),
                            'strs': come_from_del
                        }

                    # If the string comes from the insert state
                    if come_from_ins:
                        self.t_prob[I][I] = {
                            'prob': float(len(come_from_ins) /
                                          len(insert_states)),
                            'strs': list(set(come_from_ins))
                        }
                    self.t_prob[I]['strs'].extend(insert_states)

            # get emission probabilities without '.'
            num_of_dot = self.frequncy_list[i].get('.', 0)
            self.e_prob[nextM] = {
                n:
                self.frequncy_list[i][n] / (self.num_of_strings - num_of_dot)
                for n in self.frequncy_list[i] if n != '.'
            }
            i += 1

    def create_result(self):
        # write emission and transition probability
        self.t_prob = {
            n: self.t_prob[n]
            for n in self.t_prob if self.t_prob[n]['strs']
        }

        self.t_prob = sorted(self.t_prob.items(), key=_keyify)
        self.e_prob = sorted(self.e_prob.items(), key=_keyify)

        os.chdir(self.output)
        print(*self.e_prob, sep="\n",
              file=open(os.path.join(self.output, 'e.out'), 'w'))
        print(*self.t_prob, sep="\n",
              file=open(os.path.join(self.output, 't.out'), 'w'))

        #import chart
        #chart.generate_chart(self.e_prob, self.char_list)


def main():

    parser = argparse.ArgumentParser(
        description='generate profile hidden markov')
    parser.add_argument('--input', help='input file location')
    parser.add_argument('--output', help='output file location')
    args = parser.parse_args()

    model = ProfileHiddenMarkovMoldel(args.input, args.output)
    model.build_model()
    model.create_result()


if __name__ == '__main__':
    main()
