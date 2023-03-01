from multiprocessing import Process, Manager
from math import ceil


def thread_execute(func, thread_iterable, thread_num, return_dict, args):
    """
    Выполнение потока.
    :param func: функция для запуска.
    :param thread_iterable: перебираемый функцией объект, обрезанный для потока.
    :param thread_num: номер потока.
    :param return_dict: общий словарь для возвращаемых значений функций.
    :param args: аргументы функции.
    :return: возвращаемое значение функции потока.
    """
    if isinstance(args, dict):
        ret = func(thread_iterable, **args)
    elif args is None:
        ret = func(thread_iterable)
    else:
        ret = func(thread_iterable, *args)
    if return_dict is not None:
        return_dict[thread_num] = ret
    return ret


def run_function_in_multiprocess(func, iterable, args,
                                 threads: int = 2, return_type='dict', order_type='intervals'):
    """
    Запуск функции в мультипотоке. Первым аргументом в функции должен быть iterable,
    а потом остальные аргументы из args.
    :param func: функция для запуска.
    :param args: аргументы функции за исключением перебираемого объекта в виде словаря, списка или кортежа.
    :param iterable: перебираемый функцией объект.
    :param threads: число потоков.
    :param return_type: возвращаемый тип данных, если None - ничего не возвращать.
    :param order_type: тип порядка отнесения последовательности данных к процессам.
    Если intevals, то последовательность разбивается на интервалы в исходном порядке,
    иначе - в процесс N берется каждый N элемент последовательности.
    :return True - при окончании выполнения всех потоков.
    """
    if len(iterable) < threads:
        threads = len(iterable)
    step = ceil(len(iterable) / threads)
    processes = []
    manager = Manager()
    return_dict = manager.dict()
    for thread_num in range(threads):
        if order_type == 'intervals':
            if isinstance(iterable, dict):
                thread_iterable = {}
                i = -1
                for key, val in iterable.items():
                    i += 1
                    if i < thread_num * step:
                        continue
                    elif i == thread_num * step + step:
                        break
                    thread_iterable[key] = val
            else:
                thread_iterable = iterable[thread_num * step:thread_num * step + step]
        else:
            if isinstance(iterable, dict):
                thread_iterable = {}
                i = -1
                for key, val in iterable.items():
                    i += 1
                    if i == threads:
                        i = 0
                    if i == thread_num:
                        thread_iterable[key] = val
            else:
                thread_iterable = []
                i = -1
                for val in iterable:
                    i += 1
                    if i == threads:
                        i = 0
                    if i == thread_num:
                        thread_iterable.append(val)
        thread_args = (func, thread_iterable, thread_num, return_dict, args)
        processes.append(Process(target=thread_execute, args=thread_args))
        processes[thread_num].start()
    for thread_num in range(threads):
        processes[thread_num].join()
    if return_type is not None:
        if return_type in ('set', 'dict'):
            if return_type == 'set':
                return_value = set()
            else:
                return_value = {}
            for d in range(threads):
                if d in return_dict.keys():
                    val = return_dict[d]
                    return_value.update(val)
        elif return_type in ('list', 'tuple'):
            return_value = []
            for d in range(threads):
                if d in return_dict.keys():
                    val = return_dict[d]
                    return_value += val
            if return_type == 'tuple':
                return_value = tuple(return_value)
        else:
            return_value = return_dict
        return return_value
