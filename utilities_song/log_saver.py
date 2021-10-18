import datetime
import os
import glob
import shutil


class LogReader():
    """ Read in data from a night log 
    (based on the LogSaver class from Waltzcontrol)."""

    def __init__(self, current_log):
        """ Construct instance."""
        self.current_log = current_log

    def get_header_and_lines(self):
        """ Read the current log_file and return list of header lines and
            list of data lines.

            :returns: Tuple of lists of header lines and data lines.
                      Complete lines as strings with all special characters.
        """
        #try:
        with open(self.current_log, 'r') as file:
            header = [next(file) for count in range(9)]
            lines = [line for line in file]
        #except Exception:
        #    return None, None
        return header, lines

    def check_log_empty(self):
        """ Check if the current log sheet is empty.

            :returns: True if log is empty, False if data in log
        """
        _, lines = self.get_header_and_lines()

        return (lines == [])
    
    def get_saved_run(self):
        """ Get the current run saved in current log."""
        # Read in old header:
        old_header, lines = self.get_header_and_lines()

        try:
            interesting_lines = old_header[3] + old_header[4]
        except TypeError:
            return None

        index = interesting_lines.find('Run:')
        run = interesting_lines[
            index + len('Run:'):
            index + len('Run:') + 10]
        run = run.strip()

        return run

    def get_saved_date(self):
        """ Get the first and the last saved dates in the midpoint column."""
        # Read in old header:
        header, lines = self.get_header_and_lines()

        date_list = []
        try:
            for line in lines:
                date_string = line.split('|')[3]
                try:
                    date = datetime.datetime.strptime(date_string.strip(),
                                                      '%Y-%m-%dT%H:%M:%S.%f')
                    date_list.append(date)
                except Exception:
                    continue
        except TypeError:
            return None

        try:
            date_list.sort()
        except Exception:
            return None

        first = date_list[0].date()
        last = date_list[-1].date()

        result_string = self.format_header_date(first, last)

        return result_string

    def format_header_date(self, start_day, end_day):
        """ Return Header Style Date String (e.g. April 19/20 2019)

            :param start_day: Start Date as datetime.date
            :param end_day: End Date as datetime.date

            :returns: formatted date string
            :rtype: string
        """
        # If both dates have the same month: we want the format:
        # Jan 01/02 2019
        if start_day.month == end_day.month:
            result_string = "{} {}/{} {}".format(
                end_day.strftime('%b'),
                start_day.strftime('%d'),
                end_day.strftime('%d'),
                end_day.strftime('%Y'))
        # If the months are different the format should be
        # Jan 31/Feb 01 2019
        # We assume that no observation will take place at new year's eve
        # In this case the year of the end time would be saved
        else:
            result_string = "{} {}/{} {} {}".format(
                start_day.strftime('%b'),
                start_day.strftime('%d'),
                end_day.strftime('%b'),
                end_day.strftime('%d'),
                end_day.strftime('%Y'))
        return result_string

    def make_header(self, observer='', run='', telescope='',
                    date='', focus='', weather=''):
        """ Create Header as string.

            :param string observer: Name of Observer
            :param string run: Run Number
            :param string telescope: Name of Telescope
            :param string date: Date of Observation
                               (potentially with start and end date)
            :param string focus: Focus Position
            :param string weather: Weather Condition

            :returns: Header as string.

        """
        line1 = 'Waltz Spectrograph Observing Log\n'
        line2 = '\n'
        line3 = '{:-<80}\n'.format('')
        line4 = 'Observer:{:<20}  Run:{:<10}    Telescope:{:<20}\n'.format(
            observer, run, telescope)
        line5 = 'UT Date:{:<20}   Focus:{:<10}  Weather:{:<20}\n'.format(
            date, focus, weather)
        line6 = line3
        line7 = self.format_line('d*.fits', 'Object', 'I2', 'Mid-Time(UTC)',
                                 'Exp', 'Counts', 'Comment')
        line8 = self.format_line('number', 'Name', '(Y/N)', ' ',
                                 'time', ' ', ' ')
        line9 = line3

        header = (line1 + line2 + line3 + line4 + line5 +
                  line6 + line7 + line8 + line9)

        return header

    def format_line(self, d_nr='', obj_name='', i2='', mid_time='',
                    exp_time='', counts='', comment=''):
        """ Return formatted data line.

            :param string d_nr: Data Number as in d * .fits
            :param string obj_name: Name of Object e.g. (HIP12345, Wideflat)
            :param string i2: Y or N
            :param string mid_time: PMT Midpoint in UTC
            :param string exp_time: Exposure Time in seconds
            :param string counts: Maximal counts in column
            :param string comment: Comment as string

            :returns: Formatted Data Line
            :rtype: string
        """
        try:
            exp_time = str(round(float(exp_time), 3))
        except (TypeError, ValueError):
            pass

        line = """{:<8}|{:<12}|{:<5}|{:<25}|{:<9}|{:<6}|{}\n""".format(
            d_nr, obj_name, i2, mid_time, exp_time, counts, comment)

        return line

    def append_line(self, line):
        """ Append one line to the log_sheet without affecting the rest of the
            log.

            :param string line: Line to be appended to File as given
                                by format_line
        """
        try:
            with open(self.current_log, 'a') as file:
                file.write(line)
        except Exception as e:
            raise(e)
            #logger.exception("")

    def add_line(self, line):
        """ Add one line to the logfile.

            Check if data_nr exists already and replace the old line if so.
            If not append the new line.
            Split up the given line for the data_nr stored in it and search the
            current log file for this number and replace this line.
            If the data_nr does not exist yet, append it to the file.

            :param string line: Line to be appended to File as given
                                by format_line
        """
        # Get data_nr
        data_nr = line.split('|')[0]

        # Read in the whole file
        header, lines = self.get_header_and_lines()
        try:
            data_nrs = [line.split('|')[0] for line in lines]
        except TypeError:
            return None

        try:
            line_index = data_nrs.index(data_nr)
        except ValueError:
            # If the data point is not in the file yet: Just append it.
            self.append_line(line)
            return None

        # Replace the line
        lines[line_index] = line

        # And rewrite the whole file
        try:
            with open(self.current_log, 'w') as file:
                for line in header:
                    file.write(line)
                for line in lines:
                    file.write(line)
        except Exception as e:
            raise(e)
            #logger.exception("")

    def add_comment(self, data_nr, comment):
        """ Add comment to existing data line.

            :param string data_nr: Data Number that shall be edited
            :param string comment: Comment to be added
        """
        # First get old line
        _, lines = self.get_header_and_lines()
        try:
            data_nrs = [(line.split('|')[0]).strip() for line in lines]
        except TypeError:
            return None

        try:
            line_index = data_nrs.index(data_nr)
        except ValueError as e:
            raise(e)
            #logger.error('Line could not be found!')
            #return None

        old_line = lines[line_index]
        line_data = old_line.split('|')
        line_data = line_data[:6]
        line_data.append(comment)

        new_line = self.format_line(line_data[0], line_data[1],
                                    line_data[2], line_data[3],
                                    line_data[4], line_data[5],
                                    line_data[6])
        self.add_line(new_line)

    def replace_header(self, header):
        """ Replace the header of the file but leave the rest of the content.

            :param string header: formatted header as given by make_header
        """
        # Read in the whole file
        _, lines = self.get_header_and_lines()

        # And rewrite the whole file
        try:
            with open(self.current_log, 'w') as file:
                file.write(header)
                for line in lines:
                    file.write(line)
        except Exception as e:
            raise(e)
            #logger.exception("")

    def finalize_log(self, run, auto_date=True):
        """ Finalize Log: Take current log file, calculate the date and save it
            to the specified filepath.

            Filepath will be automatically determined from the run id.

            :param string run: ID of run
            :param bool auto_date: If True: Calculate the
                              date part of the header from the current date.
            :returns: New Filepath where the logsheet is stored
                      (if copying is successfull)
            :rtype: string
        """
        if auto_date:
            # We want to automatically set the date variable from the current
            # date. But we always need two dates since you observe over the
            # date transition in the night. So the idea is to check the
            # time whether it is in the morning or in the evening (since you
            # might finish before midnight). If it is in the morning, determine
            # yesterday's date and save both yesterday's dand today's date
            # If it is still in the evening, determine tomorrows date
            # and save today's and tomorrow's date
            now = datetime.datetime.now()
            if 0 < now.hour < 12:
                end_day = now.date()
                start_day = end_day - datetime.timedelta(days=1)
            else:
                start_day = now.date()
                end_day = start_day + datetime.timedelta(days=1)
            try:
                result_string = self.format_header_date(start_day, end_day)
            except Exception:
                #logger.exception("")
                result_string = ''

            self.replace_single_header_infos(date=result_string, run=run)

        # Now whether you want to automatically calculate the date or not
        # Copy the file to its new filepath
        # Parent Directory: Where all data will end up
        # For using it on Dane's PC only: (Remove that at a later stage)
        if os.uname().nodename == 'lx13':
            parent_dir = '/home/dspaeth/data/logsheets'
        else:
            parent_dir = '/home/observer/data/logsheets'

        # The file should be stored in the format: "run_id".logsheet"night_nr"
        # Where the night_nr is just the next highest number
        # that does not exist in this run
        search_format = '/{}.logsheet*'.format(run)
        existing_files = glob.glob(parent_dir + search_format)
        existing_number_files = [
            filename.split('logsheet')[-1] for filename in existing_files]
        try:
            int_existing_numbers = [int(filename)
                                    for filename in existing_number_files]
        except ValueError:
            # If you have a file that does not match the data###.fits format
            # (### indicates a number), trying to transform that to a int will
            # raise a ValueError:
            # In this case: Do not use list comprehension
            # (can't catch exception there), but use regular loop and
            # skip all names where this happens

            int_existing_numbers = []
            for filename in existing_number_files:
                try:
                    int_existing_numbers.append(
                        int(filename))
                except ValueError:
                    pass

        # Get the largest number:
        try:
            currently_largest = max(int_existing_numbers)
            next_number = currently_largest + 1
        except ValueError:
            next_number = 1

        new_file_path = '{}/{}.logsheet{}'.format(parent_dir, run, next_number)

        # And copy file
        try:
            shutil.copy(self.current_log, new_file_path)
            # Create new log to make sure that no data gets messed up
            self.new_log(auto_date=False)
            return new_file_path
        except IOError as e:
            raise(e)
            #logger.error('Unable to copy to', new_file_path)

    def replace_single_header_infos(self, observer=None, run=None,
                                    telescope=None, date=None,
                                    focus=None, weather=None):
        """ Replace one info in the header. Keep the rest of the file fixed.

            :param string observer: Name of Observer
            :param string run: Run Number
            :param string telescope: Name of Telescope
            :param string date: Date of Observation
                               (potentially with start and end date)
            :param string focus: Focus Position
            :param string weather: Weather Condition
        """

        # Read in old header:
        old_header, lines = self.get_header_and_lines()

        try:
            interesting_lines = old_header[3] + old_header[4]
        except TypeError:
            return None

        info_dict = {}
        for info, length in (('Observer', 20),
                             ('Run', 10),
                             ('Telescope', 20),
                             ('Date', 20),
                             ('Focus', 10),
                             ('Weather', 20)):
            index = interesting_lines.find(info + ':')
            variable = interesting_lines[
                index + len(info) + 1:
                index + len(info) + 1 + length]
            variable = variable.strip()
            info_dict[info] = variable

        if observer:
            info_dict['Observer'] = observer
        if run:
            info_dict['Run'] = run
        if telescope:
            info_dict['Telescope'] = telescope
        if date:
            info_dict['Date'] = date
        if focus:
            info_dict['Focus'] = focus
        if weather:
            info_dict['Weather'] = weather

        new_header = self.make_header(observer=info_dict['Observer'],
                                      run=info_dict['Run'],
                                      telescope=info_dict['Telescope'],
                                      date=info_dict['Date'],
                                      focus=info_dict['Focus'],
                                      weather=info_dict['Weather'])

        # And rewrite the whole file
        try:
            with open(self.current_log, 'w') as file:
                file.write(new_header)
                for line in lines:
                    file.write(line)
        except Exception as e:
            raise(e)
            #logger.exception("")

    
