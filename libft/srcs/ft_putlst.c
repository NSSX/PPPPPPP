/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_putlst.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/20 17:03:01 by igomez            #+#    #+#             */
/*   Updated: 2014/11/21 19:12:49 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

void	ft_putlst(t_list *lst, void (*print)(void *))
{
	if (!lst || !print)
		ft_putendl(NULL);
	else
	{
		while (lst)
		{
			print(lst->content);
			ft_putstr(" > ");
			lst = lst->next;
		}
		ft_putendl(NULL);
	}
}
